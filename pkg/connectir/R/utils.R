vcat <- function(verbose, msg, ..., newline=TRUE) {
    if (verbose) {
        cat(sprintf(msg, ...))
        if (newline) cat("\n")
    }
}

vstop <- function(msg, ...) stop(sprintf(msg, ...))

set_parallel_procs <- function(nforks=1, nthreads=1, verbose=FALSE) {
    vcat(verbose, "Setting %i parallel forks", nforks)
    suppressPackageStartupMessages(library("doMC"))
    registerDoMC()
    nprocs <- getDoParWorkers()
    if (nforks > nprocs)
        vstop("# of forks %i is greater than the actual # of processors (%i)", nforks, nprocs)
    options(cores=nforks)
    
    vcat(verbose, "Setting %i threads to be used for matrix algebra operations", nthreads)
    if (nthreads > nprocs)
        vstop("# of threads %i is greater than the actual # of processors (%i)", nthreads, nprocs)
    if (existsFunction("setMKLthreads")) {
        setMKLthreads(nthreads)
    } else {
        # cover all our bases
        Sys.setenv(MKL_NUM_THREADS=nthreads)
        Sys.setenv(GOTO_NUM_THREADS=nthreads)
        Sys.setenv(OMP_NUM_THREADS=nthreads)
    }
    
    invisible(TRUE)
}

n2gb <- function(x) x*8/1024^3
gb2n <- function(x) x/8*1024^3

# opts => list(blocksize=0, memlimit=6, verbose=TRUE)
# return opts
get_subdist_memlimit <- function(opts, nsubs, nvoxs, subs.ntpts) {
    printf <- function(msg, ..., newline=TRUE) vcat(opts$verbose, msg, ..., newline=newline)
    
    # amount of RAM used in GB for functionals
    mem_used4func <- sum(sapply(subs.ntpts, function(x) n2gb(x*nvoxs)))
    printf("...%.2f GB of memory used for functional data", mem_used4func)
    
    # amount of RAM used for distance matrix
    mem_used4dmat <- n2gb(nsubs^2 * nvoxs)*2    # *2 cuz copies are made (e.g. saving & gower)
    printf("...%.2f GB of memory used for 2 distance matrices", mem_used4dmat)
    n4onemap <- nvoxs * nsubs
    printf("...%.2f GB of memory used for one voxel's connectivity map across subjects", 
            n2gb(n4onemap))
    
    # functionals and 2 subject distances not held in memory @ same time
    mem_dat <- ifelse(mem_used4func > (mem_used4dmat/2), 
                        mem_used4func + mem_used4func, 
                        mem_used4dmat)
    
    # set blocksize if auto
    if (opts$blocksize == 0) {
        # minimum amount of RAM needed
        ## mem_used4func + memory for 2 connectivity maps per subjects
        min_mem_needed <- n2gb(n4onemap*2*getDoParWorkers()) + mem_dat
        
        # limit in RAM use
        mem_limit <- as.numeric(opts$memlimit)
        if (mem_limit < min_mem_needed) {
            vstop("You require at least %.2f GB of memory but are limited to %i GB. Please set the --memlimit option to a higher number in order to continue.", min_mem_needed, mem_limit)
        } else {
            printf("...memory limit is %.2f GB and a minimum of %.2f GB needed", 
                    mem_limit, min_mem_needed)
        }
        
        # amount of RAM for connectivity matrix
        mem_used4conn <- mem_limit - mem_used4func - mem_used4dmat

        # block size
        printf("...autosetting blocksize to -> ", newline=F)
        opts$blocksize <- floor(gb2n(mem_used4conn)/(n4onemap*getDoParWorkers()))
        printf("%i (with RAM limit of %.2f GB)", opts$blocksize, mem_limit)
    } else {
        printf("...adjusting blocksize based on # of processors and will use: ", newline=F)

        # set block size based on # of processors
        opts$blocksize <- floor(opts$blocksize/getDoParWorkers())

        # calculate amount of memory that will be used
        mem_used <- n2gb(opts$blocksize * n4onemap * getDoParWorkers()) + mem_dat
        printf("%.2f GB of RAM", mem_used)
    }
    
    return(opts)
}

bm_rowsum <- function(bigmat) {
    as.vector(.Call("bm_rowsum", bigmat, PACKAGE = "connectir"))
}

bm_rowmean <- function(bigmat) {
    as.vector(.Call("bm_rowmean", bigmat, PACKAGE = "connectir"))
}

sub.big.matrix_wrapper <- function(x, firstRow=1, lastRow=NULL, firstCol=1, lastCol=NULL, backingpath=NULL)
{
    if (is.shared(x))
        sub.big.matrix(x, firstRow, lastRow, firstCol, lastCol, backingpath)
    else
        sub.big.matrix_nonshared(x, firstRow, lastRow, firstCol, lastCol)
}

sub.big.matrix_nonshared <- function(x, firstRow=1, lastRow=NULL, firstCol=1, lastCol=NULL)
{
  rowOffset <- firstRow-1
  colOffset <- firstCol-1
  rbm <- x
  if (is.null(lastRow)) lastRow <- nrow(rbm)
  if (is.null(lastCol)) lastCol <- ncol(rbm)
  numCols <- lastCol-firstCol+1
  numRows <- lastRow-firstRow+1
  if (colOffset < 0 || rowOffset < 0 || numCols < 1 || numRows < 1 ||
      colOffset+numCols > ncol(rbm) || rowOffset+numRows > nrow(rbm))
  {
    stop(paste("A sub.big.matrix object could not be created",
               "with the specified parameters"))
  }
  .Call("SetRowOffsetInfo", rbm@address, 
        as.double(rowOffset + .Call("GetRowOffset", rbm@address)), 
        as.double(numRows), PACKAGE="bigmemory" )
  .Call("SetColumnOffsetInfo", rbm@address, 
        as.double(colOffset + .Call("GetColOffset", rbm@address)),
        as.double(numCols), PACKAGE="bigmemory")
  return(rbm)
}
