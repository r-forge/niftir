#' @nord
setGeneric('is.big.niftiXd', 
    function(x) standardGeneric('is.big.niftiXd')
)

#' @nord
setMethod('is.big.niftiXd',
    signature(x='big.niftiXd'),
    function(x) return(TRUE)
)

#' @nord
setMethod('is.big.niftiXd',
    definition=function(x) return(FALSE)
)

#' @nord
niftir.big.matrix <- function(nr, nc, output=NULL, overwrite=F, ...) {
    if (is.null(output)) {
        bf <- character()
        df <- character()
        bp <- NULL
        
        ybigmat <- big.matrix(nr, nc, ...)
    } else {
        output <- abspath(output)
        fpath <- dirname(output)
        bin_fname <- basename(output)
        desc_fname <- sprintf("%s.desc", rmext(bin_fname))
        
        bp <- fpath
        bf <- file.path(fpath, bin_fname)
        df <- file.path(fpath, desc_fname)
        
        if (file.exists(file.path(fpath, bin_fname))) {
            if (!overwrite) {
                warning(sprintf("file-backed big.matrix '%s/%s' already exists; NOT overwriting!", fpath, bin_fname))
                return(list(
                    exists = TRUE,
                    bm = attach.big.matrix(file.path(fpath, desc_fname)),
                    bp = bp,
                    bf = bf,
                    df = df
                ))
            } else {
                warning(sprintf("file-backed big.matrix '%s' already exists; WILL overwrite !", bf))
                file.remove(bf, df)
            }
        }
        
        ybigmat <- big.matrix(nr, nc, descriptorfile=desc_fname, 
                    backingfile=bin_fname, backingpath=fpath, ...)
    }
    
    return(list(
        exists = FALSE,
        bm = ybigmat,
        bp = bp,
        bf = bf,
        df = df
    ))
}

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
    signature(x="big.niftiXd", backingpath="NULL"),
    function(x, backingpath) {
        free.memory(x)
    }
)

setMethod("free.memory",
    signature(x="list", backingpath="NULL"),
    function(x, backingpath) {
        free.memory(x)
    }
)
    
#' @nord
setMethod("free.memory",
    signature(x="big.niftiXd", backingpath="missing"),
    function(x) {
        if (!is.filebacked(x))
            stop("input to free.memory cannot be a non-filebacked big.matrix")
        # free up memory
        .Call("CDestroyBigMatrix", x@address, PACKAGE="bigmemory")
        gc()
        # reattach matrix
        tmp <- attach.big.matrix(x@descriptorfile)
        x@address <- tmp@address
        # done!
        return(x)
    }
)

#' @nord
setMethod("free.memory",
    signature(x="list", backingpath="missing"),
    function(x) {
        xs <- x
        # check input
        lapply(xs, function(x) {
            if (!is.big.niftiXd(x))
                stop("input to free.memory must be bigniftiXd object")
            if (!is.filebacked(x))
                stop("input to free.memory cannot be a non-filebacked big.matrix")
        })
        # free memory
        lapply(xs, function(x) 
            .Call("CDestroyBigMatrix", x@address, PACKAGE="bigmemory")
        )
        gc()
        # reattach matrices
        for (i in 1:length(xs)) {
            tmp <- attach.big.matrix(xs[[i]]@descriptorfile)
            xs[[i]]@address <- tmp@address
        }
        # done!
        return(xs)
    }
)
