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

vbca2 <- function(bigmat, cols, ztransform=FALSE, outmat=NULL, ...) {
        A <- deepcopy(bigmat, cols=cols)
    #    A <- sub.big.matrix(bigmat, firstCol=firstCol, lastCol=lastCol)
    	B <- bigmat
    	if (is.null(outmat))
    	    C <- big.matrix(ncol(A), ncol(B), init=0, type="double", ...)
    	else if (ncol(A)==nrow(outmat) && ncol(B)==ncol(outmat))
    	    C <- outmat
    	else
    	    stop("dimensions of outmat are out of whack")
    	.Call("big_cor", A, B, C, PACKAGE = "connectir")
    	if (ztransform)
    	     .Call("big_ztransform", C, PACKAGE = "connectir")
    	invisible(C)
}

vbca_batch2 <- function(subs.bigmats, cols, ztransform=FALSE, ...) {
    nsubjects <- length(subs.bigmats)
    lapply(subs.bigmats, function(x) vbca2(x, cols, ztransform, ...))
}
