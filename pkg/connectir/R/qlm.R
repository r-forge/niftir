# Quick Multiple Linear Regression using Big Matrices

# 1. Have a set of big matrices, for one voxel @ a time calculate a linear regression

# Given a design matrix => check it's rank; calculate the dd first + save
# want function that converts a bm => arma mat

setOldClass(c("qfit", "qcontrasts"))

formula_to_mat <- function(formula, data) {
    mf <- model.frame(formula=formula, data=data)
    X <- model.matrix(attr(mf, "terms"), data=mf)
    return(X)
}

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

# Residuals
setGeneric('qlm_residuals', 
    function(y, X, check.rank=FALSE, residuals=NULL, ...)
        standardGeneric('qlm_residuals')
)

setMethod('qlm_residuals',
    signature(y='big.matrix', X='big.matrix'),
    function(y, X, check.rank=FALSE, residuals=NULL, ...) {
        n = nrow(X); k = ncol(X); m = ncol(y)
        
        if (check.rank) {
            Xrank <- qlm_rank(X)
            if (Xrank < k)
                stop("design matrix is rank deficient.")
        }
        
        if (is.null(residuals)) {
            residuals = big.matrix(nrow(y), ncol(y), type="double", ...)
        } else if (nrow(residuals) != n && ncol(residuals) != m) {
            stop("size mismatch for output residuals")
        }
        
        .Call("big_qlm_residuals", y, X, residuals, PACKAGE = "connectir")
        
        return(residuals)
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

# Checks design matrix
# and if rank deficient removes offending columns
check_design_matrix <- function(X, to.quit=FALSE) {
    qrhs <- qr(X)
    if (qrhs$rank < ncol(X)) {
        if (to.quit)
            stop("model is rank deficient")
        
        vcat(verbose, " warning: model is rank deficient")
        bad_cols <- paste(qrhs$pivot[-c(1:qrhs$rank)], collapse=", ")
        bad_colnames <- paste(colnames(X)[bad_cols], collapse=", ")
        vcat(verbose, " will drop cols %s (%s) from design matrix", 
                        bad_cols, bad_colnames)
        
        X <- rhs[, X$pivot, drop = FALSE]
        X <- rhs[, 1:X$rank, drop = FALSE]

    }
    
    return(X)
}

vox_glm <- function(funclist, evs, cons, blocksize, outmats, bp=NULL, 
                    verbose=TRUE, parallel=FALSE, shared=parallel, 
                    ztransform=FALSE)
{
    vcat(verbose, "...setup")
    nsubs <- length(funclist)
    nvoxs <- ncol(funclist[[1]])
    voxs <- 1:nvoxs
    ncons <- nrow(cons)
    blocks <- niftir.split.indices(1, nvoxs, by=blocksize)
    progress <- ifelse(verbose, "text", "none")
    k <- qlm_rank(evs)
    if (k < ncol(evs))
        stop("EV matrix is rank deficient")
    dd <- qlm_dd(evs)
    if (!is.shared(outmats[[1]]))
        stop("outmats must be shared")
    if (is.filebacked(outmats[[1]]) && is.null(bp))
        stop("backingpath (bp) must be given if outmats are file-baced")
    if (length(outmats) != ncons)
        stop("cons doesn't match with outmats")
    
    gfun <- function(bi) {
        vcat(verbose, "...block %i", bi)
        
        first <- blocks$starts[bi]; last <- blocks$ends[bi]
        sub_voxs <- first:last
        sub_nvoxs <- length(sub_voxs)
        
        tmp_outmats <- lapply(1:ncons, function(x) {
                            big.matrix(nvoxs, sub_nvoxs, init=0, type="double", shared=shared)
                        })
        
        # correlate
        vcat(verbose, "....correlate")
        subs.cormaps <- vbca_batch2(funclist, c(first, last), 
                                    ztransform=ztransform, 
                                    type="double", shared=shared)
        
        # loop through each seed
        vcat(verbose, '....glm')
        sfun <- function(si) {
            # combine and transpose correlation maps for 1 seed
            seedCorMaps <- big.matrix(nsubs, nvoxs-1, type="double", shared=FALSE)
            .Call("subdist_combine_and_trans_submaps", subs.cormaps, as.double(si), 
                  as.double(voxs[-sub_voxs[si]]), seedCorMaps, PACKAGE="connectir")
                        
            # glm
            fit <- qlm_fit(seedCorMaps, evs, shared=FALSE)
            res <- qlm_contrasts(fit, cons, dd, shared=FALSE)
            
            for (ci in 1:ncons)
                tmp_outmats[[ci]][-voxs[sub_voxs[si]],si] <- res$tvalues[ci,]
            
            rm(seedCorMaps, fit, res); gc(FALSE)
            
            return(NULL)
        }
        
        llply(1:sub_nvoxs, sfun, .parallel=parallel, .progress=progress, .inform=T)
        
        # copy over temp data
        vcat(verbose, '....save')
        for (ci in 1:ncons) {
            sub_outmat <- sub.big.matrix(outmats[[ci]], firstCol=first, lastCol=last)
            deepcopy(x=tmp_outmats[[ci]], y=sub_outmat)
            if (!is.null(bp)) {
                flush(sub_outmat); flush(outmats[[ci]])
                outmats[[ci]] <- free.memory(outmats[[ci]], bp)
            }
        }
        
        rm(subs.cormaps, tmp_outmats); gc(FALSE)
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    for (bi in 1:blocks$n)
        gfun(bi)
    
    invisible(outmats)
}
