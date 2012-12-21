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

# copied from bigmemory
becleanupcols <- function(cols=NULL, nc=NULL, colnames=NULL) {
  if (is.null(cols)) cols <- 1:nc
  else {
    if (!is.numeric(cols) & !is.character(cols) & !is.logical(cols))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(cols))
      if (is.null(colnames)) stop("column names do not exist.")
      else cols <- mmap(cols, colnames)
    if (is.logical(cols)) {
      if (length(cols) != nc)
        stop(paste("column vector length must match the number of",
                   "columns of the matrix."))
      cols <- which(cols)
    }
    tempj <- .Call("CCleanIndices", as.double(cols), as.double(nc), PACKAGE="bigmemory")
    if (is.null(tempj[[1]])) stop("Illegal column index usage in extraction.\n")
    if (tempj[[1]]) cols <- tempj[[2]]
  }
  return(cols)
}

# copied from bigmemory
becleanuprows <- function(rows=NULL, nr=NULL, rownames=NULL) {
  if (is.null(rows)) rows <- 1:nr
  else {
    if (!is.numeric(rows) & !is.character(rows) & !is.logical(rows))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(rows))
      if (is.null(rownames)) stop("row names do not exist.")
      else rows <- mmap(rows, rownames)
    if (is.logical(rows)) {
      if (length(rows) != nr)
        stop(paste("row vector length must match the number of",
                   "rows of the matrix."))
      rows <- which(rows)
    }
    tempj <- .Call("CCleanIndices", as.double(rows), as.double(nr), PACKAGE="bigmemory")
    if (is.null(tempj[[1]])) stop("Illegal row index usage in extraction.\n")
    if (tempj[[1]]) rows <- tempj[[2]]
  }
  return(rows)
}

bedeepcopy <- function(x, x.cols=NULL, x.rows=NULL, 
                       y=NULL, y.cols=NULL, y.rows=NULL, 
                       type=NULL, separated=NULL, 
                       backingfile=NULL, backingpath=NULL,
                       descriptorfile=NULL, shared=TRUE)
{
    x.cols <- becleanupcols(x.cols, ncol(x), colnames(x))
    x.rows <- becleanuprows(x.rows, nrow(x), rownames(x))
    if (nrow(x) > 2^31-1)
      stop(paste("Too many rows to copy at this point in time;",
                 "this may be fixed in the future."))
    if (is.null(type)) type <- typeof(x)
    if (is.big.matrix(x)) {
      if (is.null(separated)) separated <- is.separated(x)
    } else {
      separated <- FALSE
    }
    if (is.null(y)) {
      y <- big.matrix(nrow=length(x.rows), ncol=length(x.cols), type=type, init=NULL,
                    dimnames=dimnames(x), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile, shared)
    }
    if (typeof(x) != type)
        stop("The type of x and the type of y must be the same")
    if (is.separated(x) != separated)
        stop("x and y must have both have separated columns or both have not separated columns")
    y.cols <- becleanupcols(y.cols, ncol(y), colnames(y))
    y.rows <- becleanuprows(y.rows, nrow(y), rownames(y))
    if (is.big.matrix(x) && is.big.matrix(y)) {
        .Call("BigDeepCopyMain", x@address, as.double(x.rows), as.double(x.cols), 
              y@address, as.double(y.rows), as.double(y.cols), 
              PACKAGE="bigextensions")
    } else {
        for (i in 1:length(cols)) y[,i] <- x[rows,cols[i]]
    }
    return(y)
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
    signature(x="list", backingpath="NULL"),
    function(x, backingpath) {
        free.memory(x)
    }
)

