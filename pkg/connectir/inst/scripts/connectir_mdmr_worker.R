suppressPackageStartupMessages(library("optparse"))


# Make option list
option_list <- list(
    make_option(c("-i", "--indir"), type="character", default=NULL, help="Input subdist directory (required)", metavar="subdist"),
    make_option(c("-f", "--formula"), type="character", default=NULL, help="a typical R model formula that specifies the factors or continuous variables that may expain the variance in each voxel's subject distance matrix", metavar="'A + B:C'"),
    make_option(c("-m", "--model"), type="character", default=NULL, help="Filename of a comma separated file with participant info in R friendly format where column names correspond to formula values... (required)", metavar="csv"),
    make_option("--whichsubs", type="character", default=NULL, help="Filename with a list of subject indices to use from the subject distance matrices (default is to use all of them)", metavar="text-file"),
    make_option("--expr", type="character", default=NULL, help="An expression based on the model that is used to restrict the subjects examined (can either use this or --whichsubs, not both)", metavar="expression"),
    make_option("--strata", type="character", default=NULL, help="Only compute permutations within groups, you can specify the name of a column in your '--model' that indicates these groups (optional)", metavar="name"),
    make_option(c("-p", "--permutations"), type="integer", default=4999, help="Number of permutations to conduct for each voxel [default: %default]", metavar="number"),
    make_option("--factors2perm", type="character", default=NULL, help="Which factors (e.g., A and B) to permute from the formula specified [default: all of them]", metavar="'A,B'"),
    make_option(c("-c", "--forks"), type="integer", default=2, help="Number of computer processors to use in parallel by forking the complete processing stream [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel by multi-threading matrix algebra operations [default: %default]", metavar="number"),
    make_option("--blocksize", type="integer", default=0, dest="blocksize", help="How many sets of voxels should used in each iteration of computing the pseudo F-statistics (0 = auto) [default: %default]", metavar="number"),
    make_option("--memlimit", type="integer", default=6, dest="memlimit", help="If blocksize is set to auto (--blocksize=0), then will set the blocksize to use a maximum of RAM specified by this option  [default: %default]", metavar="number"),
    make_option("--save-perms", action="store_true", default=FALSE, dest="saveperms", help="Save all the permuted psuedo-F stats? [default: %default]"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output if it already exists (default is not to overwrite already existing output)"),
    make_option(c("-d", "--debug"), action="store_true", default=TRUE, help="Like verbose but will also print more helpful error messages when --forks is >1"), 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output")
)
    
# Make class/usage
parser <- OptionParser(usage = "%prog [options] output", option_list=option_list, add_help_option=TRUE)

# Parse
parser_out <- parse_args(parser, positional_arguments = TRUE)
args <- parser_out$args
opts <- parser_out$options

# Check options/arguments
if (length(args) != 1) {
    print_help(parser)
    quit(save="no", status=1)
}

saved_opts <- list(args=args, opts=opts)

set_parallel_procs <- function(nforks=1, nthreads=1, verbose=FALSE) {
    vcat(verbose, "Setting %i parallel forks", nforks)
    suppressPackageStartupMessages(library("blasctl"))
    suppressPackageStartupMessages(library("doMC"))
    registerDoMC()
    nprocs <- getDoParWorkers()
    if (nforks > nprocs) {
        vstop("# of forks %i is greater than the actual # of processors (%i)", 
              nforks, nprocs)
    }
    options(cores=nforks)
    
    vcat(verbose, "Setting %i threads for matrix algebra operations", 
         nthreads)
    nprocs <- omp_get_max_threads()
    if (nthreads > nprocs) {
        vstop("# of threads %i is greater than the actual # of processors (%i)", 
              nthreads, nprocs)
    }
    
    if (existsFunction("setMKLthreads")) {
        vcat(verbose, "...using Intel's MKL")
        setMKLthreads(nthreads)
    } else {
        # cover all our blases
        vcat(verbose, "...using GOTOBLAS or Other")
        blas_set_num_threads(nthreads)
        omp_set_num_threads(nthreads)
    }
    
    # Not sure if these env vars matter?
    Sys.setenv(OMP_NUM_THREADS=nthreads)
    Sys.setenv(GOTO_NUM_THREADS=nthreads)
    Sys.setenv(MKL_NUM_THREADS=nthreads)
    Sys.setenv(OMP_DYNAMIC=TRUE)
    Sys.setenv(MKL_DYNAMIC=TRUE)
    
    invisible(TRUE)
}

tryCatch({
  
  # load connectir
  suppressWarnings(suppressPackageStartupMessages(library("connectir")))

  # parallel processing setup
  set_parallel_procs(opts$forks, opts$threads, opts$verbose)  
  # use foreach parallelization and shared memory?
  parallel_forks <- ifelse(opts$forks == 1, FALSE, TRUE)
  
  ###
  # Check Inputs
  ###
  vcat(opts$verbose, "Checking options")
  
  # check variables
  if (is.null(opts$indir))
      stop("must specify -i/--indir option")
  if (is.null(opts$model))
      stop("must specify -m/--model option")
  if (is.null(opts$formula))
      stop("must specify -f/--formula option")
  if (opts$overwrite)
      stop("Right now the overwrite function isn't implemented")
  if (opts$permutations < 2)
      stop("# of permutations must be greater than 1")
  if (!is.null(opts$expr) && !is.null(opts$whichsubs))
      stop("cannot specify both --expr and --whichsubs")
  
  # check paths exist
  opts$indir <- abspath(opts$indir)
  invisible(check_subdist(opts$indir))
  if (!file.exists(opts$model))
      stop("-m/--model ", opts$model, " does not exist")
  opts$outdir <- file.path(opts$indir, args[1])
  if (file.exists(opts$outdir))
      stop("Output mdmr directory '", opts$outdir,  "' already exists")
  if (!file.exists(dirname(opts$outdir)))
      stop("Basepath to mdmr directory '", dirname(opts$outdir), "' doesn't exist")
  
  
  ###
  # Setup inputs
  ###
  vcat(opts$verbose, "Setting up inputs")
  
  if (opts$debug) {
      verbosity <- 2
  } else if (opts$verbose) {
      verbosity <- 1
  } else {
      verbosity <- 0
  }
  
  # data dimensions
  vcat(opts$verbose, "...data dimensions")
  tmp <- attach.big.matrix(file.path(opts$indir, "subdist.desc"))
  nsubs <- sqrt(nrow(tmp))
  nvoxs <- ncol(tmp)
  rm(tmp); gc(FALSE, TRUE)

  # formula
  vcat(opts$verbose, "...formula")
  formula <- as.formula(sprintf(". ~ %s", opts$formula))
  vars <- as.character(as.list(attr(terms(formula), "variables"))[-c(1:2)])
  
  # factors2perm
  vcat(opts$verbose, "...factors to permute")
  if (!is.null(opts$factors2perm)) {
      vcat(opts$verbose, "...factors2perm")
      
      opts$factors2perm <- sub(", ", ",", opts$factors2perm)
      opts$factors2perm <- strsplit(opts$factors2perm, ",")[[1]]
      
      for (x in opts$factors2perm) {
          if (!(x %in% vars))
              vstop("Factor to permute '%s' not found in formula", x)
      }
      
      nfactors <- length(opts$factors2perm)
  } else {
      nfactors <- length(attr(terms(formula), "term.labels"))
  }
  
  # model
  vcat(opts$verbose, "...model")
  model <- read.csv(opts$model)
  ## checks
  for (v in vars) {
      if (is.null(model[[v]]))
          vstop("Factor '%s' doesn't match any column in model file", x)
  }
  
  # strata
  if (!is.null(opts$strata)) {
      vcat(opts$verbose, "...strata")
      if (is.null(model[[opts$strata]]))
          stop("Strata given but doesn't match any column in the model file")
      else
          opts$strata <- model[[opts$strata]]
  }
  
  # filter subjects
  if (!is.null(opts$whichsubs)) {
      vcat(opts$verbose, "...whichsubs")
      filter.subs <- TRUE
      which.subs <- as.numeric(read.table(opts$whichsubs)[,1])
      if (all(which.subs==0 || which.subs==1))
          which.subs <- which(which.subs==1)
      if (length(which.subs) == 0)
          stop("no subjects left to analyze based on --whichsubs")
  } else if (!is.null(opts$expr)) {
      vcat(opts$verbose, "...expr")
      filter.subs <- TRUE
      which.subs <- eval(parse(text=sprintf("with(model, which(%s))", opts$expr)))
      if (length(which.subs) == 0)
          stop("no subjects left to analyze based on --expr")
  } else {
      filter.subs <- FALSE
  }
  ## model
  if (filter.subs) {
      if (nrow(model) == nsubs) {
          model <- model[which.subs,]
      } else if (nrow(model) != length(which.subs)) {
          stop(paste("# of rows in model don't match # of subjects in", 
                     "distance matrix or # of filtered subjects"))
      }
  } else {
      if (nrow(model) != nsubs)
          stop("# of rows in model file don't match # of subjects in distance matrix")
  }
  
  # output
  vcat(opts$verbose, "...creating output directory '%s'", opts$outdir)
  dir.create(opts$outdir)
  
  # subject distances
  if (filter.subs) {
      fname <- file.path(opts$indir, "subdist.desc")
      xdist <- filter_subdist_fb(fname, which.subs, opts$outdir, opts$memlimit, 
                                 parallel=parallel_forks, verbose=opts$verbose, 
                                 gower=TRUE)
      nsubs <- length(which.subs)
      xdist.path <- opts$outdir
  } else {
      vcat(opts$verbose, "Reading in subject distances")
      fname <- file.path(opts$indir, "subdist_gower.desc")
      xdist <- attach.big.matrix(fname)
      xdist.path <- opts$indir
  }
  
  # check
  vcat(opts$verbose, "...checking input")
  tmp <- matrix(xdist[,1], nsubs, nsubs)
  check_gmat(tmp)
  rm(tmp); invisible(gc(FALSE, TRUE))
  
  
  ###
  # Memory Demands: superblocksize & blocksize
  ###
  
  nperms <- opts$permutations
  opts <- get_mdmr_memlimit(opts, nsubs, nvoxs, nperms, nfactors)
  
  
  ###
  # Compute MDMR
  ###
  start.time <- Sys.time()
  
  if (opts$saveperms) {
      fperms.path <- opts$outdir
  } else {
      fperms.path <- NULL
  }
  
  res.mdmr <- mdmr(xdist, formula, model, nperms=opts$permutations, 
                   superblocksize=opts$superblocksize, blocksize=opts$blocksize, 
                   strata=opts$strata, factors2perm=opts$factors2perm, 
                   verbose=verbosity, parallel=parallel_forks, 
                   G.path=xdist.path, fperms.path=fperms.path)
  rm(xdist)
  invisible(gc(FALSE, TRUE))
  
  end.time <- Sys.time()
  vcat(opts$verbose, "MDMR is done! It took: %.2f minutes\n", 
       as.numeric(end.time-start.time, units="mins"))
  
  
  ###
  # Save MDMR Results
  ###
  save_mdmr(res.mdmr, opts$indir, opts$outdir, formula, opts$verbose)
  rm(res.mdmr)
  invisible(gc(FALSE, TRUE))

}, warning = function(ex) {
  cat("\nA warning was detected: \n")
  cat(ex$message, "\n\n")
  cat("Called by: \n")
  print(ex$call)
  cat("\nSaving options...\n")
  save(saved_opts$args, saved_opts$opts, file="called_options.rda")
}, error = function(ex) {
  cat("\nAn error was detected: \n")
  cat(ex$message, "\n\n")
  cat("Called by: \n")
  print(ex$call)
  cat("\nSaving options...\n")
  save(saved_opts$args, saved_opts$opts, file="called_options.rda")
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
