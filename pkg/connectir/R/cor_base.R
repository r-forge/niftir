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
