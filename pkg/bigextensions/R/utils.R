deepcopy <- function (inmat, cols = NULL, rows = NULL, outmat=NULL, ...) {
    if (!is.big.matrix(inmat)) 
        stop("input must be a big.matrix")
    if (is.null(cols)) 
        cols <- 1:ncol(inmat)
    if (is.null(rows)) 
        rows <- 1:nrow(inmat)
    if (is.null(outmat)) {
        # note: type and seperated cannot be user supplied
        outmat <- big.matrix(length(rows), length(cols), type = typeof(inmat), 
            separated = is.separated(inmat), ...)
    }
    
    .Call("BigDeepCopyMain", inmat@address, outmat@address, as.double(rows), 
        as.double(cols))
    return(outmat)
}

#' @nord
setMethod('diag', 
    signature(x="big.matrix"),
    function(x) {
        len <- min(dim(x))
        .Call("GetDiagMain", x@address, as.double(len))
    }
)

setMethod('diag<-', 
    signature(x="big.matrix"),
    function(x, value) {
        if (length(value) == 1)
            values <- rep(value, min(dim(x)))
        else if (length(value) == min(dim(x)))
            values <- value
        else
            stop("replacement diagonal has wrong length")
        .Call("SetDiagMain", x@address, as.double(values))
        return(x)
    }
)

#' @nord
setGeneric("free.memory",
    function(x, backingpath, ...)
        standardGeneric('free.memory')
)

setMethod("free.memory",
    signature(x="big.matrix", backingpath="character"),
    function(x, backingpath) {
        if (!is.filebacked(x))
            return(x)
        # free up memory
        d <- describe(x)
        .Call("CDestroyBigMatrix", x@address, PACKAGE="bigmemory")
        gc()
        # reattach matrix
        tmp <- attach.big.matrix(d, backingpath=backingpath)
        x@address <- tmp@address
        # done!
        return(x)
    }
)

setMethod("free.memory",
    signature(x="list", backingpath="NULL"),
    function(x, backingpath) {
        free.memory(x)
    }
)

