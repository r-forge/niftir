vbca <- function(bigmat, cols, ztransform=FALSE, outmat=NULL, ...) {
    A <- deepcopy(bigmat, cols=cols)
#    A <- sub.big.matrix(bigmat, firstCol=firstCol, lastCol=lastCol)
	B <- bigmat
	if (is.null(outmat))
	    C <- big.matrix(ncol(A), ncol(B), init=0, type="double", ...)
	else if (ncol(A)==nrow(outmat) && ncol(B)==ncol(outmat))
	    C <- outmat
	else
	    stop("dimensions of outmat are out of whack")
	ALPHA <- 1/(nrow(B) - 1)
	dgemm(C=C, A=A, B=B, TRANSA='t', ALPHA=ALPHA)
	if (ztransform)
	    atanh(C)
	invisible(C)
}

vbca_batch <- function(subs.bigmats, cols, ztransform=FALSE, ...) {
    nsubjects <- length(subs.bigmats)
    lapply(1:nsubjects, function(i) vbca(subs.bigmats[[i]], cols, ztransform, ...))
}


## NEW CODE USING ARMADILLO

big_cor <- function(x, y=NULL, z=NULL, byrows=FALSE, 
                    x_firstCol=1, x_lastCol=ncol(x), 
                    y_firstCol=NULL, y_lastCol=NULL, 
                    z_firstCol=1, z_lastCol=NULL, 
                    ...) 
{
    # Setup 'y'
    if (is.null(y)) {
        y <- x
        if (is.null(y_firstCol))
            y_firstCol <- x_firstCol
        if (is.null(y_lastCol))
            y_lastCol <- x_lastCol
    }
    if (is.null(y_firstCol))
        y_firstCol <- 1
    if (is.null(y_lastCol))
        y_lastCol <- ncol(y)
    
    # Setup/Check 'z'
    if (is.null(z)) {   # Create output matrix
        if (byrows) {
            nrow <- nrow(x)
            ncol <- nrow(y)
        } else {
            nrow <- x_lastCol - x_firstCol + 1
            ncol <- y_lastCol - y_firstCol + 1
        }
        z <- big.matrix(nrow, ncol, ...)
        if (z_firstCol!=1 || !is.null(z_lastCol))
            warning("user specifications of z_firstCol and z_lastCol not used (must give z)", 
                    immediate.=TRUE)
        z_firstCol <- 1; z_lastCol <- ncol
    } else if (is.null(z_lastCol)) {
        z_lastCol <- ncol(z)
    }
    
    # Check all inputs
    if (!is.big.matrix(x) || !is.big.matrix(y) || !is.big.matrix(z))
        stop("inputs must be big matrices")
    
    # How to perform the crossproduct?
    if (byrows) {
        method <- "big_tcor"
    } else {
        method <- "big_cor"
    }
    
    # Correlate!
    .Call(method, x, y, z, 
            as.double(x_firstCol), as.double(x_lastCol), 
            as.double(y_firstCol), as.double(y_lastCol), 
            as.double(z_firstCol), as.double(z_lastCol), 
            PACKAGE = "connectir")
    
    invisible(z)
}

big_ztransform <- function(x) {
    if (!is.big.matrix(x))
        stop("input must be a big matrix")
    .Call("big_ztransform", x, PACKAGE="connectir")
}

# incols => c(firstCol, lastCol)
# outcols => c(firstCol, lastCol)
vbca2 <- function(inmat, incols=c(1,ncol(inmat)), ztransform=FALSE, 
                  outmat=NULL, outcols=c(1,incols[2]), ...) 
{
	if (length(incols)!=2 || length(outcols)!=2)
	    stop("incols and outcols must be a vector of length 2: c(firstCol, lastCol)")
	
	# TODO: figure out better way to handle this
	if (is.null(outmat))
	    z_lastCol = NULL
	else
	    z_lastCol = outcols[2]
	
	outmat <- big_cor(x=inmat, z=outmat, 
	                  x_firstCol=incols[1], x_lastCol=incols[2], 
	                  y_firstCol=1, y_lastCol=ncol(inmat), 
	                  z_firstCol=outcols[1], z_lastCol=z_lastCol)
	if (ztransform) atanh(outmat)
	
	invisible(outmat)
}

vbca_batch2 <- function(subs.bigmats, cols, ztransform=FALSE, ...) {
    nsubjects <- length(subs.bigmats)
    lapply(subs.bigmats, function(x) vbca2(x, cols, ztransform, ...))
}

