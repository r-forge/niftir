#!/usr/bin/env Rscript

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
    make_option(c("-i", "--indir"), type="character", default=NULL, help="Input subdist directory (required)", metavar="subdist"),
    make_option(c("-o", "--outdir"), type="character", default="mdmr", help="Directory to output MDMR results (this will be within the input subdist directory)", metavar="MDMR"),
    make_option(c("-m", "--model"), type="character", default=NULL, help="Filename of a comma separated file with participant info in R friendly format where column names correspond to formula values... (required)", metavar="csv"),
    make_option("--usesubs", type="character", default=NULL, help="Filename with a list of subject indices to use from the subject distance matrices (default is to use all of them)", metavar="text-file"),
    make_option("--strata", type="character", default=NULL, help="Only compute permutations within groups, you can specify the name of a column in your '--model' that indicates these groups (optional)", metavar="name"),
    make_option(c("-p", "--permutations"), type="integer", default=4999, help="Number of permutations to conduct for each voxel [default: %default]", metavar="number"),
    make_option("--factors2perm", type="character", default=NULL, help="Which factors to permute from the formula specified [default: all of them]", metavar="comma-separated list"),
    make_option(c("-c", "--cores"), type="integer", default=2, help="Number of computer processors to use in parallel [default: %default]", metavar="number"),
    make_option(c("-t", "--threads"), type="integer", default=1, help="Number of computer processors to use in parallel for MKL library [default: %default]", metavar="number"),
    make_option("--blocksize", type="integer", default=0, dest="blocksize", help="How many sets of voxels should used in each iteration of computing the pseudo F-statistics (0 = auto) [default: %default]", metavar="number"),
    make_option("--overwrite", action="store_true", default=FALSE, help="Overwrite output if it already exists (default is not to overwrite already existing output)"),
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE, help="Print extra output [default]"),
    make_option(c("-q", "--quiet"), action="store_false", dest="verbose", help="Print little output")
)
    
# Make class/usage
parser <- OptionParser(usage = "%prog [options] formula", option_list=option_list, add_help_option=TRUE)

# Parse
parser_out <- parse_args(parser, positional_arguments = TRUE)
args <- parser_out$args
opts <- parser_out$options

# Check options/arguments
if (length(args) < 1) {
    print_help(parser)
    quit(save="no", status=1)
}

###
# Parallel processing setup
###
printf("03a. Setting %i cores to be used", opts$cores)
if (opts$cores > 1) {
    printf("...setting parallel processing with doMC")
    suppressPackageStartupMessages(library("doMC"))
    registerDoMC()
    if (opts$cores > getDoParWorkers())
    	stop("Number of -c/--cores specified '", opts$cores, "' is greater than the actual number of cores '", getDoParWorkers(), "'")
}
options(cores=opts$cores)
#if (existsFunction("setMKLthreads")) {
	printf("03b. Setting %i MKL threads to be used", opts$threads)
	printf("...setting number of threads for MKL")
	Sys.setenv(MKL_NUM_THREADS=opts$threads)
#}

###
# Check Inputs
###
printf("01. Checking options")
## check variables
if (is.null(opts$indir))
    stop("must specify -i/--indir option")
if (is.null(opts$model))
    stop("must specify -m/--model option")
if (opts$overwrite)
    stop("Right now the overwrite function isn't implemented")
if (opts$permutations < 2)
    stop("# of permutations must be greater than 1")
## load connectir
suppressWarnings(suppressPackageStartupMessages(library("connectir")))
## check paths exist
opts$indir <- abspath(opts$indir)
invisible(check_subdist(opts$indir))
if (!file.exists(opts$model))
    stop("-m/--model ", opts$model, " does not exist")
opts$outdir <- file.path(opts$indir, opts$outdir)
if (file.exists(opts$outdir))
    stop("Output mdmr directory '", opts$outdir,  "' already exists")
if (!file.exists(dirname(opts$outdir)))
    stop("Basepath to mdmr directory '", dirname(opts$outdir), "' doesn't exist")

###
# Setup inputs
###
printf("02. Setting up inputs")

# formula
printf("...formula")
formula <- as.formula(sprintf(". ~ %s", paste(args, collapse=" ")))

# factors2perm
if (!is.null(opts$factors2perm)) {
    printf("...factors2perm")
    opts$factors2perm <- sub(", ", ",", opts$factors2perm)
    opts$factors2perm <- strsplit(opts$factors2perm, ",")[[1]]
}

# model
printf("...model")
model <- read.csv(opts$model)
if (!is.null(opts$strata)) {
    if (is.null(model[[opts$strata]]))
        stop("Strata given but doesn't match any column in the input model")
    else
        opts$strata <- model[[opts$strata]]
}
# TODO: check if below works with more than one permuted factor
if (!is.null(opts$factors2perm) && is.null(model[[opts$factors2perm]]))
    stop("Factors to permute given but doesn't match any column in the input model")

# subject distances (gower's centered)
printf("...subject distances")
## read
tmp <- attach.big.matrix(file.path(opts$indir, "subdist_gower.desc"))
## copy into memory
if (is.null(opts$usesubs))
    opts$usesubs <- 1:sqrt(nrow(tmp))
else
    opts$usesubs <- as.numeric(read.table(opts$usesubs)[,1])
xdist <- slice.subdist(tmp, subs=opts$usesubs)  # will create a copy
rm(tmp); invisible(gc(FALSE))
## check
tmp <- matrix(xdist[,1], sqrt(nrow(xdist)), sqrt(nrow(xdist)))
check_gmat(tmp)
rm(tmp); invisible(gc(FALSE))


###
# Block size
###
printf("04. Dealing with the block size and RAM")

# functions
n2gb <- function(x) x*8/1024^3
gb2n <- function(x) x/8*1024^3

# general info
nsubs <- sqrt(nrow(xdist))
nvoxs <- ncol(xdist)

# RAM for distance matrix
mem_used4dmat <- n2gb(nsubs^2 * nvoxs)

# RAM for Fperms matrix
mem_used4fperms <- opts$permutations * nvoxs

if (opts$blocksize==0) {
    printf("...autosetting blocksize to -> ", newline=F)
    
    # minimum amount of RAM needed (assume blocksize of 2)
    min_mem_needed <- n2gb(
        opts$permutations * 2 * 2   # for temporary matrices
        + nsubs^2 * 2
    ) * getDoParWorkers() + mem_used4dmat + mem_used4fperms
    
    # limit in RAM use
    mem_limit <- Sys.getenv("CONNECTIR_RAM_LIMIT")
    mem_limit <- ifelse(mem_limit == "", 6, as.numeric(mem_limit))
    if (mem_limit < min_mem_needed)
        stop(sprintf("You require at least %.2f GB of memory but are limited to %i GB. Please set the environmental variable CONNECTIR_RAM_LIMIT to a higher number to continue.", min_mem_needed, mem_limit))
    
    # amount of RAM for mdmr
    mem_used4mdmr <- mem_limit - mem_used4dmat - mem_used4fperms
    
    # block size
    opts$blocksize <- floor(gb2n(
        mem_used4mdmr/((2*opts$permutations+nsubs^2)*getDoParWorkers())
    ))
    printf("%i", opts$blocksize)
    
    # clear variables
    rm(min_mem_needed, mem_limit, mem_used4mdmr)
    
}

# temporary RAM needed while computing MDMR
## error variance and explained variance
mem_used4temp <- n2gb(opts$permutations * opts$blocksize * 2 * getDoParWorkers())

# copy of subdist
mem_used4chunk <- n2gb(nsubs^2 * opts$blocksize * getDoParWorkers())

# total ram
mem_used <- mem_used4dmat + mem_used4fperms + mem_used4temp + mem_used4chunk

printf("...will be using %.2f GB of RAM", mem_used)
rm(nsubs, nvoxs, mem_used4dmat, mem_used4fpers, mem_used4temp, mem_used4chunk, mem_used)
invisible(gc(FALSE))


###
# Compute MDMR
###
printf("05. Computing MDMR")
start.time <- Sys.time()

res.mdmr <- mdmr(xdist, formula, model, nperms=opts$permutations, factors2perm=opts$factors2perm, block.size=opts$blocksize, verbose=opts$verbose, strata=opts$strata)

end.time <- Sys.time()
printf("MDMR is done! It took: %.2f minutes\n", as.numeric(end.time-start.time, units="mins"))


###
# Save MDMR Results
###
printf("06. Saving MDMR Results")
save_mdmr(res.mdmr, opts$indir, opts$outdir, formula, verbose=opts$verbose)
