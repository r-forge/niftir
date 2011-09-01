suppressPackageStartupMessages(library("optparse"))

# General Function(s)
printf <- function(msg, ..., newline=TRUE) {
    if (opts$verbose) {
        cat(sprintf(msg, ...))
        if (newline) cat("\n")
    }
}

# Make option list
option_list <- list(
    make_option(c("-i", "--infuncs"), type="character", default=NULL, dest="infuncs", help="File containing paths of different 4D functional images in one column (required)", metavar="file"),
    make_option(c("-m", "--inmasks"), type="character", default=NULL, dest="inmasks", help="File containing paths of different 3D masks for each functional image (required). Make sure that this list of masks are in the same order as the -i/--infuncs file.", metavar="file"),
    make_option("--ztransform", action="store_true", default=FALSE, dest="ztransform", help="Fischer Z-Transform the correlations before calculating the distance between participants"),
    make_option("--seedmask", type="character", default=NULL, help="Mask to select the voxels that will be used to correlate with each voxel in the rest of the brain (or anything within the specified --brainmask)", metavar="file"),
    make_option("--brainmask", type="character", default=NULL, help="When computing each whole-brain connectivity map, this mask will restrict which parts of the whole-brain are to be considered", metavar="file"),
    make_option("--blocksize", type="integer", default=0, help="How many sets of voxels should be used in each iteration of computing the correlation (0 = auto) [default: %default]", metavar="number"),
    make_option("--memlimit", type="integer", default=6, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option(c("-c", "--cores"), type="integer", default=1, help="Number of computer processors to use in parallel [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel for MKL library [default: %default]", metavar="number"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output that already exists (default is not to overwrite already existing output)"),
    make_option("--no-link-functionals", action="store_true", default=FALSE, help="Will not create soft links to each of the functional images with the subdist directory"),
    make_option("--check-functionals", action="store_true", default=FALSE, help="This will check if any of the input functional images has a point with a value equal to or less than zero, which is not allowed"),
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output")
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

  ###
  # Parallel processing setup
  ###
  printf("04. Setting %i cores to be used", opts$cores)
  if (opts$cores > 1) {
      printf("...setting parallel processing with doMC")
      suppressPackageStartupMessages(library("doMC"))
      registerDoMC()
      if (opts$cores > getDoParWorkers())
      	stop("Number of -c/--cores specified '", opts$cores, "' is greater than the actual number of cores '", getDoParWorkers(), "'")
  }
  options(cores=opts$cores)
	
	printf("04. Setting %i MKL threads to be used", opts$threads)
	printf("...setting number of threads for MKL")
  if (existsFunction("setMKLthreads")) {
  	setMKLthreads(opts$threads)
  } else {
  	Sys.setenv(MKL_NUM_THREADS=opts$threads)
  }
  
  # want shared memory or not
  use_shared = ifelse(opts$cores == 1, FALSE, TRUE)
  
  suppressWarnings(suppressPackageStartupMessages(library("connectir")))

  start.time <- Sys.time()

  ###
  # Check Arguments
  ###
  printf("01. Checking required inputs")
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
  ## get input functionals
  infiles <- sapply(as.character(read.table(opts$infuncs)[,1]), function(fp) {
      if (!file.exists(fp))
          stop("One of the input functionals does not exist: ", fp)
      abspath(fp)
  })
  n <- length(infiles)
  ## get input masks
  inmasks <- sapply(as.character(read.table(opts$inmasks)[,1]), function(fp) {
      if (!file.exists(fp))
          stop("One of the input functionals does not exist: ", fp)
      abspath(fp)
  })
  if (length(inmasks) != n)
      stop("Number of masks is not the same as the number of functional images")

  ###
  # Check Options
  ###
  printf("02. Checking optional inputs")
  if (!is.null(opts$seedmask)) {
      if(!file.exists(opts$seedmask))
          stop("--seedmask file ", opts$seedmask, " does not exist")
      opts$seedmask <- abspath(opts$seedmask)
  }
  if (!is.null(opts$brainmask)) {
      if(!file.exists(opts$brainmask))
          stop("--brainmask file ", opts$rowmask, " does not exist")
      opts$brainmask <- abspath(opts$brainmask)
  }


  ###
  # Read in inputs
  ###
  printf("05. Setting up inputs")
  ## remove existing output
  if (opts$overwrite)
      stop("Right now the overwrite function isn't implemented")

  ## masks
  if (is.null(opts$brainmask)) {
      prebrainmask <- NULL
  } else {
      printf("...reading brain mask")
      prebrainmask <- read.mask(opts$brainmask)
  }
  if (is.null(opts$seedmask)) {
  	if (!is.null(prebrainmask))
  		preseedmask <- prebrainmask
  	else
      	preseedmask <- NULL
  } else {
      printf("...reading seed mask")
      preseedmask <- read.mask(opts$seedmask)
  }

  ## overlap mask
  printf("...creating overlap of masks across participants")
  maskoverlap <- create_maskoverlap(inmasks)
  
  ## seed mask
  printf("...creating final seed mask")
  if (!is.null(preseedmask)) {
      if (length(preseedmask) != length(maskoverlap))
          stop("length of seedmask and maskoverlap not the same")
  #    if (sum(preseedmask[!maskoverlap]) > 0)
  #        warning(sprintf("Seed mask '%s' contains some voxels that don't overlap across all participants", opts$seedmask))
      seedmask <- preseedmask & maskoverlap
  } else {
      seedmask <- maskoverlap
  }

  ## brainmask
  printf("...creating final brain mask")
  if (!is.null(prebrainmask)) {
      if (length(prebrainmask) != length(maskoverlap))
          stop("length of brainmask and maskoverlap not the same")
  #    if (sum(prebrainmask[!maskoverlap]) > 0)
  #        warning(sprintf("Brain mask '%s' contains some voxels that don't overlap across all participants", opts$brainmask))
      brainmask <- prebrainmask & maskoverlap
  } else {
      brainmask <- maskoverlap
  }
  if (!all(seedmask[brainmask]==TRUE))
      stop("For now the brainmask must contain all elements of the seedmask")
  
  # set block size and memlimit
  nsubs <- length(infiles)
  nvoxs <- sum(brainmask)
  ntpts <- sapply(infiles, function(x) {
      hdr <- read.nifti.header(x)
      if (length(hdr$dim) != 4)
        stop("Input functional file must be 4 dimensions: ", x, ", but is ", length(hdr$dim))
      return(hdr$dim[[4]])
  })
  opts <- get_memlimit(opts, nsubs, nvoxs, ntpts)
  
  ## functional data
  printf("...reading and masking the functional data")
  # 1. check if this is ok
  funclist <- load_and_mask_func_data(infiles, brainmask, shared=use_shared, type="double", 
                                      verbose=opts$verbose)
  invisible(gc(FALSE))
  if (opts$"check-functionals") {
    printf("...checking functional data")
    for (i in 1:length(funclist)) {
      func <- funclist[[i]]
      x <- colmin(func)
      if (any(x <= 0)) {
          stop("Masked functional data has points with a value equal or less than zero, this is not allowed: ", infiles[[i]])
      }
    }
  }

  # create the subdist directory (and get the subject distance matrix)
  printf("...creating subdist directory and files")
  masks <- list(maskoverlap=maskoverlap, seedmask=seedmask, brainmask=brainmask)
  if (!is.null(preseedmask)) masks$preseedmask <- preseedmask
  if (!is.null(prebrainmask)) masks$prebrainmask <- prebrainmask
  # note: this subdist isn't filebacked yet
  subdist <- create_subdist(outdir, infiles, masks, opts) 
  
  # clear memory of unneeded stuff + get seedinds
  seedinds <- which(seedmask[brainmask])
  rm(masks, maskoverlap, preseedmask, prebrainmask, seedmask)   # save brainmask
  invisible(gc(F))

  end.time <- Sys.time()
  printf("Setup is done. It took: %.2f minutes\n", as.numeric(end.time-start.time, units="mins"))


  ###
  # Compute the subdist
  ###
  start.time <- Sys.time()
  
  printf("06. Computing subject distances")
  compute_subdist(funclist, subdist, seed_inds=seedinds, blocksize=opts$blocksize, ztransform=opts$ztransform, verbose=opts$verbose)
  rm(funclist)
  invisible(gc(FALSE))

  end.time <- Sys.time()
  printf("Distance computation is done! It took: %.2f minutes\n", as.numeric(end.time-start.time, units="mins"))

  ###
  # Save the subdist
  ###
  printf("07. Saving subject distances")
  tmp <- deepcopy(subdist, backingpath=outdir, backingfile="subdist.bin", descriptorfile="subdist.desc")
  rm(tmp)
  invisible(gc(FALSE))

  if (any(is.na(subdist[2,]))) {
      print(subdist[2,])
      stop("Found NA's in second row of subdist")
  }
  
  printf("...saving 3D image zcheck.nii.gz")
  zcheck <- (subdist[2,]!=0)*1 + 1
  hdr <- read.nifti.header(infiles[1])
  hdr$dim <- hdr$dim[1:3]; hdr$pixdim <- hdr$pixdim[1:3]
  write.nifti(zcheck, hdr, brainmask, outfile=file.path(outdir, "zcheck.nii.gz"))
  rm(brainmask)
  if (any(zcheck==1))
    fprintf("There are some bad voxels...see zcheck.nii.gz")
  
  
  ###
  # Create gower matrix
  ###
  start.time <- Sys.time()

  printf("08. Creating gower's centered matrices")
  gdist <- gower.subdist(subdist)
  rm(subdist)
  invisible(gc(FALSE))
  printf("...checking")
  mat <- matrix(gdist[,1], sqrt(nrow(gdist)))
  check_gmat(mat)
  
  end.time <- Sys.time()
  printf("Centering of matrices done! It took: %.2f minutes\n", as.numeric(end.time-start.time, units="mins"))


  ###
  # Save gower matrix
  ###
  printf("09. Saving gower's centered matrices")
  tmp <- deepcopy(gdist, backingpath=outdir, backingfile="subdist_gower.bin", descriptorfile="subdist_gower.desc")
  rm(tmp, gdist)
  invisible(gc(FALSE))

}, warning = function(ex) {
  cat("\nA warning was detected: \n")
  cat(ex$message, "\n\n")
  cat("Called by: \n")
  print(ex$call)
}, error = function(ex) {
  cat("\nAn error was detected: \n")
  cat(ex$message, "\n\n")
  cat("Called by: \n")
  print(ex$call)
  cat("\nSaving options...\n")
  save(args, opts, printf, file="called_options.rda")
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
