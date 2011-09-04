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
    
    # set blocksize if auto
    if (opts$blocksize == 0) {
        # minimum amount of RAM needed
        ## mem_used4func + memory for 2 connectivity maps per subjects
        min_mem_needed <- n2gb(n4onemap*2*getDoParWorkers()) + mem_used4func + mem_used4dmat
        
        # limit in RAM use
        mem_limit <- as.numeric(opts$memlimit)
        if (mem_limit < min_mem_needed)
            stop(sprintf("You require at least %.2f GB of memory but are limited to %i GB. Please set the --memlimit option to a higher number in order to continue.", min_mem_needed, mem_limit))
        else
            printf("...memory limit is %.2f GB", mem_limit)
        
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
        mem_used <- n2gb(opts$blocksize * n4onemap)
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
