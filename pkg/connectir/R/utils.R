vcat <- function(verbose, ...) {
    cat(...)
    cat("\n")
}

n2gb <- function(x) x*8/1024^3
gb2n <- function(x) x/8*1024^3

# opts => list(blocksize=0, memlimit=6, verbose=TRUE)
# return opts
get_memlimit <- function(opts, nsubs, nvoxs, subs.ntpts) {
    printf <- function(msg, ..., newline=TRUE) {
        if (opts$verbose) {
            cat(sprintf(msg, ...))
            if (newline) cat("\n")
        }
    }

    # amount of RAM used in GB for functionals
    mem_used4func <- sum(sapply(subs.ntpts, function(x) n2gb(x*nvoxs)))
    
    # amount of RAM used for distance matrix
    mem_used4dmat <- n2gb(nsubs^2 * nvoxs)
    n4onemap <- nvoxs * nsubs
    
    # set blocksize if auto
    if (opts$blocksize == 0) {
        printf("...autosetting blocksize to -> ", newline=F)

        # minimum amount of RAM needed
        ## mem_used4func + memory for 2 connectivity maps per subjects
        min_mem_needed <- n2gb(n4onemap*2*getDoParWorkers()) + mem_used4func + mem_used4dmat
        
        # limit in RAM use
        mem_limit <- as.numeric(opts$memlimit)
        if (mem_limit < min_mem_needed)
            stop(sprintf("You require at least %.2f GB of memory but are limited to %i GB. Please set the --memlimit option to a higher number in order to continue.", min_mem_needed, mem_limit))
        
        # amount of RAM for connectivity matrix
        mem_used4conn <- mem_limit - mem_used4func - mem_used4dmat

        # block size
        opts$blocksize <- floor(gb2n(mem_used4conn)/(n4onemap*getDoParWorkers()))
        printf("%i (with RAM limit of %.2f GB)", opts$blocksize, mem_limit)

        # clear variables
        rm(min_mem_needed, mem_used4conn, mem_limit)
    } else {
        printf("...adjusting blocksize based on # of processors and will use: ", newline=F)

        # set block size based on # of processors
        opts$blocksize <- floor(opts$blocksize/getDoParWorkers())

        # calculate amount of memory that will be used
        mem_used <- n2gb(opts$blocksize * n4onemap)
        printf("%.2f GB of RAM", mem_used)
        rm(mem_used)
    }
    rm(n2gb, gb2n, mem_used4func, mem_used4dmat, n4onemap)
    
    return(opts)
}


