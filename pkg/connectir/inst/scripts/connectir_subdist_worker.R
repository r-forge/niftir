suppressPackageStartupMessages(library("optparse"))

# Make option list
option_list <- list(
    make_option(c("-i", "--infuncs"), type="character", default=NULL, dest="infuncs", help="File containing paths of different 4D functional images in one column (required)", metavar="file"),
    make_option(c("-m", "--inmasks"), type="character", default=NULL, dest="inmasks", help="File containing paths of different 3D masks for each functional image (required). Make sure that this list of masks are in the same order as the -i/--infuncs file.", metavar="file"),
    make_option("--ztransform", action="store_true", default=FALSE, dest="ztransform", help="Fischer Z-Transform the correlations before calculating the distance between participants"),
    make_option("--brainmask", type="character", default=NULL, help="When computing each whole-brain connectivity map, this mask will restrict which parts of the whole-brain are to be considered", metavar="file"),
    make_option("--bg", type="character", default=NULL, help="Background image (e.g., MNI152 standard brain) upon which later results might be overlaid", metavar="file"), 
    make_option("--blocksize", type="integer", default=0, help="How many sets of voxels should be used in each iteration of computing the correlation (0 = auto) [default: %default]", metavar="number"),
    make_option("--memlimit", type="integer", default=6, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option("--forks", type="integer", default=1, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option("--threads", type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output that already exists (default is not to overwrite already existing output)"),
    make_option("--no-link-functionals", action="store_true", default=FALSE, help="Will not create soft links to each of the functional images with the subdist directory"),
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output"), 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-d", "--debug"), action="store_true", default=TRUE, help="Like verbose but will also print more helpful error messages when --forks is >1")
)

# Make class/usage
parser <- OptionParser(usage = "%prog [options] output-directory", option_list=option_list, add_help_option=TRUE)

# Parse
parser_out <- parse_args(parser, positional_arguments = TRUE)
args <- parser_out$args
opts <- parser_out$options

if (length(args) < 1) {
    print_help(parser)
    quit(save="no", status=1)
}

tryCatch({

  # parallel processing setup
  set_parallel_procs(opts$forks, opts$threads, opts$verbose)  
  # use foreach parallelization and shared memory?
  parallel_forks = ifelse(opts$cores == 1, FALSE, TRUE)
  
  # load connectir
  suppressWarnings(suppressPackageStartupMessages(library("connectir")))


  ###
  # Check/Setup Required Inputs
  ###
  
  start.time <- Sys.time()
  
  # Check required
  vcat(opts$verbose, "Checking required inputs")
  outdir <- abspath(args[1])
  if (file.exists(outdir) && !opts$overwrite)
      stop("Output directory '", outdir, "' already exists, you can use --overwrite")
  if (is.null(opts$infuncs))
      stop("You must specify the -i/--infuncs option")
  if (is.null(opts$inmasks))
      stop("You must specify the -m/--inmasks option")
  if (!file.exists(opts$infuncs))
      stop("The file specified by -i/--infuncs does not exist")
  if (!file.exists(opts$inmasks))
      stop("The file specified by -m/--inmasks does not exist")
  
  # Prepare input functional filenames
  infiles <- sapply(as.character(read.table(opts$infuncs)[,1]), function(fp) {
      if (!file.exists(fp))
          stop("One of the input functionals does not exist: ", fp)
      abspath(fp)
  })
  n <- length(infiles)
  
  # Prepare input mask filenames
  inmasks <- sapply(as.character(read.table(opts$inmasks)[,1]), function(fp) {
      if (!file.exists(fp))
          stop("One of the input functionals does not exist: ", fp)
      abspath(fp)
  })
  if (length(inmasks) != n)
      stop("Number of masks is not the same as the number of functional images")
  
  
  ###
  # Check/Setup Optional Inputs
  ###

  vcat(opts$verbose, "Checking optional inputs")
  for (optname in c("brainmask", "bg")) {
      arg <- opts[[optname]]
      if(!file.exists(arg))
          vstop("--%s file '%s' does not exist", optname, arg)
      opts[[optname]] <- abspath(arg)
  }
  
  if (opts$debug) {
      verbosity <- 2
  } else if (opts$verbose) {
      verbosity <- 1
  } else {
      verbosity <- 0
  }
  
  
  ###
  # Read/Setup Masks
  ###
  
  vcat(opts$verbose, "Setting up masks")
  
  ## remove existing output
  if (opts$overwrite)
      stop("Right now the overwrite function isn't implemented")
  
  ## brainmask 1
  if (is.null(opts$brainmask)) {
      prebrainmask <- NULL
  } else {
      vcat(opts$verbose, "...reading brain mask")
      prebrainmask <- read.mask(opts$brainmask)
  }
  
  ## overlap mask
  vcat(opts$verbose, "...creating overlap of masks across participants")
  maskoverlap <- create_maskoverlap(inmasks)
  
  ## brainmask 2
  vcat(opts$verbose, "...creating final brain mask")
  if (is.null(prebrainmask)) {
      brainmask <- maskoverlap
  } else {
      if (length(prebrainmask) != length(maskoverlap))
          stop("length of brainmask and maskoverlap not the same")
      brainmask <- prebrainmask & maskoverlap
  }
  
  
  ###
  # Set Memory Demands
  ###
  
  nsubs <- length(infiles)
  nvoxs <- sum(brainmask)
  ntpts <- sapply(infiles, function(x) {
      hdr <- read.nifti.header(x)
      if (length(hdr$dim) != 4) {
          vstop("Input functional file '%s' must be 4 dimensions but is %i dimensional", 
                x, length(hdr$dim))
      }
      return(hdr$dim[[4]])
  })
  opts <- get_subdist_memlimit(opts, nsubs, nvoxs, ntpts)
  
  
  ###
  # Read/Prepare Functional Data
  ###
  
  vcat(opts$verbose, "Loading and masking functional data")
  reader <- gen_big_reader("nifti4d")
  load_and_mask_func_data2(infiles, reader, mask=brainmask, 
                           verbose=opts$verbose, parallel=parallel_forks, 
                           type="double", shared=parallel_forks)
  invisible(gc(FALSE, TRUE))
  
  
  ###
  # Creating output directory
  ###
  
  vcat(opts$verbose, "Creating output directory and files")
  masks <- list(maskoverlap=maskoverlap, brainmask=brainmask)
  if (!is.null(prebrainmask)) masks$prebrainmask <- prebrainmask
  dists_list <- create_subdist(outdir, infiles, masks, opts, shared=parallel_forks)
  
  # clear memory of unneeded stuff + get seedinds
  rm(masks, maskoverlap)   # save brainmask
  invisible(gc(FALSE, TRUE))
  
  end.time <- Sys.time()
  vcat(opts$verbose, "Setup is done. It took: %.2f minutes\n", 
       as.numeric(end.time-start.time, units="mins"))
  
  
  ###
  # Compute the subdists
  ###
  start.time <- Sys.time()
  
  vcat(opts$verbose, "Computing subject distances")
  checks <- compute_subdist_wrapper(funclist, dists_list, 
                                    opts$blocksize, opts$superblocksize, 
                                    verbose=verbosity, parallel=parallel_forks, 
                                    ztransform=opts$ztransform, method="pearson")  
    
  vcat(opts$verbose, "...saving zchecks")  
  hdr <- read.nifti.header(infiles[1])
  hdr$dim <- hdr$dim[1:3]; hdr$pixdim <- hdr$pixdim[1:3]
  write.nifti(checks$sdist, hdr, brainmask, odt="char", 
              outfile=file.path(outdir, "zcheck_subdist.nii.gz"))
  write.nifti(checks$gdist, hdr, brainmask, odt="char", 
              outfile=file.path(outdir, "zcheck_subdist_gower.nii.gz"))

  end.time <- Sys.time()
  vcat(opts$verbose, "Done! Total computation time: %.1f minutes\n", 
       as.numeric(end.time-start.time, units="mins"))  
  
}, warning = function(ex) {
  cat("\nA warning was detected: \n")
  cat(ex$message, "\n\n")
  cat("Called by: \n")
  print(ex$call)
  cat("\nSaving options...\n")
  save(args, opts, file="called_options.rda")
}, error = function(ex) {
  cat("\nAn error was detected: \n")
  cat(ex$message, "\n\n")
  cat("Called by: \n")
  print(ex$call)
  cat("\nSaving options...\n")
  save(args, opts, file="called_options.rda")
}, interrupt = function(ex) {
  cat("\nKill signal sent. Trying to clean up...\n")
  rm(list(ls()))
  gc(FALSE)
  cat("...success\n")
}, finally = {
  cat("\nRemoving everything from memory\n")
  rm(list=ls())
  gc(FALSE)
  cat("...sucesss\n")
})
