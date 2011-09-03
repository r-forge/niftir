# Quick Multiple Linear Regression using Big Matrices

# 1. Have a set of big matrices, for one voxel @ a time calculate a linear regression

# Given a design matrix => check it's rank; calculate the dd first + save
# want function that converts a bm => arma mat

setOldClass(c("qfit", "qcontrasts"))

# Check if design matrix is rank deficient
setGeneric('qlm_rank', 
    function(X)
        standardGeneric('qlm_rank')
)

setMethod('qlm_rank',
    signature(X='big.matrix'),
    function(X) {
        .Call("big_qlm_rank", X, PACKAGE = "connectir")
    }
)


# Get diag( solve(t(X) %*% X) )
setGeneric('qlm_dd', 
    function(X)
        standardGeneric('qlm_dd')
)

setMethod('qlm_dd',
    signature(X='big.matrix'),
    function(X) {
        as.vector(.Call("big_qlm_dd", X, PACKAGE = "connectir"))
    }
)

# Solves for y ~ X
setGeneric('qlm_fit', 
    function(y, X, check.rank=TRUE, ...)
        standardGeneric('qlm_fit')
)

setMethod('qlm_fit',
    signature(y='big.matrix', X='big.matrix'),
    function(y, X, check.rank=FALSE, outputs=NULL, ...) {
        n = nrow(X); k = ncol(X); m = ncol(y)
        
        if (check.rank) {
            Xrank <- qlm_rank(X)
            if (Xrank < k)
                stop("design matrix is rank deficient.")
        }
        
        if (is.null(outputs)) {
            outputs <- list(
                coefficients = big.matrix(ncol(X), ncol(y), type="double", ...), 
                residuals = big.matrix(nrow(y), ncol(y), type="double", ...), 
                mse = big.matrix(ncol(y), 1, type="double", ...)
            )
        } else {
            if (!all(c("coefficients", "residuals", "mse") %in% outputs))
                stop("outputs must be list of 'coefficients', 'residuals', and 'mse'.")
            if (nrow(outputs$coef) != n && ncol(outputs$coef) != m)
                stop("size mismatch for output coefficients")
            if (nrow(outputs$res) != n && ncol(outputs$res) != m)
                stop("size mismatch for output residuals")
            if (nrow(outputs$mse) != m && ncol(outputs$mse) != 1)
                stop("size mismatch for output mse")
        }
        
        .Call("big_qlm_fit", y, X, outputs$coef, outputs$res, outputs$mse, 
                as.double(n), as.double(k), as.double(m), PACKAGE = "connectir")
        
        outputs$n <- n; outputs$k <- k; outputs$m <- m
        outputs$call <- match.call()
        class(outputs) <- "qfit"
        
        return(outputs)
    }
)

# Gets tvals for given contrasts from qlm_fit results
setGeneric('qlm_contrasts', 
    function(inputs, contrasts, dd, ...)
        standardGeneric('qlm_contrasts')
)

setMethod('qlm_contrasts',
    signature(inputs='qfit', contrasts='matrix', dd='vector'),
    function(inputs, contrasts, dd, outputs=NULL, ...) {
        if (ncol(contrasts) != inputs$k)
            stop("column mismatch for contrasts")
        if (length(dd) != inputs$k)
            stop("length mismatch for dd")
        n <- nrow(contrasts)
        
        if (is.null(outputs)) {
            outputs <- list(
                coefficients = big.matrix(n, inputs$m, type="double", ...), 
                standard_errors = big.matrix(n, inputs$m, type="double", ...), 
                tvalues = big.matrix(n, inputs$m, type="double", ...)
            )
        } else {
            if (!all(c("coefficients", "standard_errors", "tvalues") %in% outputs))
                stop("outputs must be list of 'coefficients', 'standard_errors', and 'tvalues'")
            if (nrow(outputs$coef) != n && ncol(outputs$coef) != inputs$m)
                stop("size mismatch for output coefficients")
            if (nrow(outputs$st) != n && ncol(outputs$st) != inputs$m)
                stop("size mismatch for output standard errors")
            if (nrow(outputs$tval) != n && ncol(outputs$tval) != inputs$m)
                stop("size mismatch for output t-values")
        }
        
        .Call("big_qlm_contrasts", inputs$coef, inputs$mse, contrasts, dd, outputs$coef, 
                outputs$st, outputs$tval, as.double(inputs$m), PACKAGE = "connectir")
        
        outputs$n <- n
        outputs$call <- match.call()
        class(outputs) <- "qcontrasts"
        
        return(outputs)
    }
)
