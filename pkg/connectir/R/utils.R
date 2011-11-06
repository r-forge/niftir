vcat <- function(verbose, msg, ..., newline=TRUE) {
    if (verbose) {
        cat(sprintf(msg, ...))
        if (newline) cat("\n")
    }
}

vstop <- function(msg, ...) stop(sprintf(msg, ...))

vsystem <- function(verbose, cmd, ...) {
    vcat(verbose, cmd, ...)
    ret <- system(sprintf(msg, ...))
    if (ret != 0)
        stop("command failed")
}

set_parallel_procs <- function(nforks=1, nthreads=1, verbose=FALSE) {
    vcat(verbose, "Setting %i parallel forks", nforks)
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
    #nprocs <- omp_get_max_threads()
    if (nthreads > nprocs) {
        vstop("# of threads %i is greater than the actual # of processors (%i)", 
              nthreads, nprocs)
    }
    
    if (existsFunction("setMKLthreads", where=topenv(.GlobalEnv))) {
        vcat(verbose, "...using Intel's MKL")
        setMKLthreads(nthreads)
    } else {
        # cover all our blases
        vcat(verbose, "...using GOTOBLAS or Other")
        suppressPackageStartupMessages(library("blasctl"))
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

big_rowsum <- function(bigmat, firstCol=1, lastCol=ncol(bigmat)) {
    as.vector(.Call("big_rowsum", bigmat, as.double(firstCol), as.double(lastCol), 
              PACKAGE="connectir"))
}

big_rowmean <- function(bigmat, firstCol=1, lastCol=ncol(bigmat)) {
    as.vector(.Call("big_rowmean", bigmat, as.double(firstCol), as.double(lastCol), 
              PACKAGE="connectir"))
}

# x = a*x + b
big_add_multiply_scalar <- function(x, a=1, b=0, firstCol=1, lastCol=ncol(x)) {
    .Call("big_add_multiply_scalar", x, as.double(a), as.double(b), 
          as.double(firstCol), as.double(lastCol), PACKAGE="connectir")
    invisible(x)
}
