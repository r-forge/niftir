load_func_data <- function(fnames) lapply(fnames, read.big.nifti4d)

create_maskoverlap_fromfuncs <- function(funcs) {
    allmasks <- sapply(funcs, function(x) colmin(x) != 0)
    apply(allmasks, 1, all)
}

mask_func_data <- function(funcs, mask) {
    dat <- lapply(funcs, function(x) {
        x <- do.mask(x, mask)
        y <- scale(x, to.copy=TRUE) # copying to remove association w/ big.nifti4d
        rm(x)
        y
    })
    gc(FALSE)
    return(dat)
}

create_maskoverlap <- function(mask_fnames) {
    overlap <- read.mask(mask_fnames[[1]])
    for (i in 1:length(mask_fnames))
        overlap <- overlap & read.mask(mask_fnames[[i]])
    return(overlap)
}

# check, will check that dims are all same
load_and_mask_func_data <- function(fnames, mask, check=TRUE, type=NULL, verbose=FALSE, ...) {
    if (is.character(mask))
        mask <- read.mask(mask)
    if (is.character(fnames) && !is.vector(fnames))
        fnames <- list(fnames)
    
    if (opts$verbose)
        progress="text"
    else
        progress="none"
    
    dat.list <- llply(fnames, function(f) {
        if (is.null(type))
            x <- read.big.nifti4d(f, ...)
        else
            x <- read.big.nifti4d(f, type=type, ...)
        x <- do.mask(x, mask)
        y <- scale(x, to.copy=T, ...) # fix with update of deepcopy!
        rm(x)
        gc(FALSE)
        y
    }, .progress=progress)
    
    if (check) {
        nr <- nrow(dat.list[[1]])
        nc <- ncol(dat.list[[1]])
        for (i in 1:length(dat.list)) {
            #if (nr != nrow(dat.list[[i]]))
            #    warning("loading functional data...not all inputs have the same # of 'timepoints'")
            if (nc != ncol(dat.list[[i]]))
                stop("loading functional data...not all inputs have the same # of voxels")
        }
    }
    
    gc(FALSE)
    
    if (length(fnames) == 1)
        return(dat.list[[1]])
    else
        return(dat.list)
}

# BORROWED CODE FROM MASS PACKAGE
ginv <- function (X, tol = sqrt(.Machine$double.eps)) {
    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X)) 
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X)) 
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive)) 
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive)) 
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
        t(Xsvd$u[, Positive, drop = FALSE]))
}
# END BORROWED CODE

# BELOW: CODE DIRECTLY TAKEN FROM PACKAGE vegan
permuted.index <- function (n, strata) 
{
    if (missing(strata) || is.null(strata)) 
        out <- sample(n, n)
    else {
        out <- 1:n
        inds <- names(table(strata))
        for (is in inds) {
            gr <- out[strata == is]
            if (length(gr) > 1) 
                out[gr] <- sample(gr, length(gr))
        }
    }
    out
}
# ABOVE: CODE DIRECTLY TAKEN FROM PACKAGE vegan

