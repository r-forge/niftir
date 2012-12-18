#' Generate the hat matrices based on the model/formula
#' as well as associated data like degrees of freedom
#' 
#' @author Zarrar Shehzad
#' @param formula A typical model formula such as Y ~ A + B*C but where Y
#'                isn't actually used
#' @param model Data frame containing your factors (e.g., A, B, & C)
#' @param contr.unordered, contr.ordered contrasts used for the design matrix
#' @param factors2perm factors to permute
#' @param verbose boolean
#' @return list of RHS, QR of RHS, full-model H, H2s, H1s, IH, and dfs
mdmr_model <- function(formula, model, 
                        contr.unordered="contr.sum", 
                        contr.ordered="contr.poly", 
                        factors2perm=NULL, 
                        verbose=FALSE)
{
    vcat(verbose, "Generating hat matrices and like (with mdmr_model)")
    
    vcat(verbose, "...creating right-hand design matrix")
    rhs <- mdmr_model.rhs(formula, model, contr.unordered, contr.ordered)
    
    vcat(verbose, "...calculating QR decomposition")
    qrhs <- mdmr_model.qr(rhs, factors2perm)
    
    vcat(verbose, "...checking/adjusting for rank deficiencies")
    rhs <- mdmr_model.rank(rhs, qrhs)
    
    vcat(verbose, "...creating the hat matrices (H2s and IHs)")
    hats <- mdmr_model.hat_matrices(rhs, qrhs)
    
    vcat(verbose, "...calculating degrees of freedom")
    dfs <- mdmr_model.df(qrhs)
    
    vcat(verbose, "...combining elements to return")
    ret <- list(
        formula = formula, 
        model = model, 
        rhs = rhs, 
        qrhs = qrhs, 
        H2s = hats$H2s, 
        IHs = hats$IHs, 
        df.res = dfs$res, 
        df.exp = dfs$exp
    )
    
    return(ret)
}

#' INTERNAL: Create right-hand/design matrix
#' 
#' @author Zarrar Shehzad
#' @param formula A typical model formula such as Y ~ A + B*C but where Y
#'                isn't actually used;
#' @param model Data frame containing your factors (e.g., A, B, & C)
#' @param contr.unordered, contr.ordered contrasts used for the design matrix
#' @return matrix
mdmr_model.rhs <- function(formula, model, contr.unordered="contr.sum",  
                           contr.ordered="contr.poly") 
{
    n <- nrow(model)
    
    rhs.frame <- model.frame(formula, model, drop.unused.levels = TRUE)
    if (nrow(rhs.frame) != nrow(model))
        vstop("One of your factors has an empty element")
    op.c <- options()$contrasts
    options(contrasts = c(contr.unordered, contr.ordered))
    rhs <- model.matrix(formula, rhs.frame)
    options(contrasts = op.c)
    
    # Add factor.names to right-hand matrix
    factor.names <- attr(attr(rhs.frame, "terms"), "term.labels")
    attr(rhs, "factor.names") <- factor.names
    
    return(rhs)
}

#' INTERNAL: Checks that right-hand matrix isn't rank-deficient 
#' by getting the QR decomposition of the RHS
#' 
#' @author Zarrar Shehzad
#' @param rhs right hand matrix
#' @param factors2perm factors to permute
#' @return QR decomposition of rhs, this includes added attributes: 
#'         grps, u.grps, factor.names, & factor2perm
mdmr_model.qr <- function(rhs, factors2perm=NULL)
{
    qrhs <- qr(rhs)
    
    # Relate each regressor to factor id
    grps <- attr(rhs, "assign")
    grps <- grps[qrhs$pivot][1:qrhs$rank]
    u.grps <- unique(grps)
    attr(qrhs, "grps") <- grps
    attr(qrhs, "u.grps") <- u.grps
    
    # Factor names
    factor.names <- attr(rhs, "factor.names")[u.grps]
    attr(qrhs, "factor.names") <- factor.names
    
    # Permuted factors
    cnames <- colnames(rhs)
    if (is.null(factors2perm)) {
        factors2perm <- 1:length(cnames)
        factors2perm <- factors2perm[cnames!="(Intercept)"]
    }
    if (is.character(factors2perm)) {
        factors2perm <- sapply(factors2perm, function(x) 
                                which(cnames==x))
    }
    factor2perm.names <- cnames[factors2perm]
    names(factors2perm) <- factor2perm.names
    attr(qrhs, "factors2perm") <- factors2perm
    
    return(qrhs)
}

#' INTERNAL: Checks that right-handed model matrix isn't rank-deficient 
#' and adjusts the model matrix ordering and number of columns
#' 
#' @author Zarrar Shehzad
#' @param rhs Right hand matrix
#' @param qrhs QR decomposition of right hand matrix
#' @param throw_error Whether or not to stop the program if rank-deficient
#' @return QR decomposition of rhs, this includes added attributes: 
#'         grps, u.grps, factor.names, & factor2perm
mdmr_model.rank <- function(rhs, qrhs, throw_error=TRUE)
{
    # Check rank
    if (qrhs$rank < ncol(rhs)) {
        bad_cols <- paste(qrhs$pivot[-c(1:qrhs$rank)], collapse=", ")
        bad_colnames <- paste(colnames(rhs)[bad_cols], collapse=", ")
        msg <- sprintf("model is rank deficient, check cols %s (%s)", 
                        bad_cols, bad_colnames)
        if (throw_error) {
            stop(msg)
        } else {
            message(msg)
        }
    }
    
    # Fix RHS based on rank and pivot
    rhs <- rhs[, qrhs$pivot, drop = FALSE]
    rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
    
    return(rhs)
}

##' INTERNAL: Calculates the full-model hat matrix
##' 
##' @author Zarrar Shehzad
##' @param rhs right-hand matrix
##' @return matrix
#mdmr_model.hat_matrix_full <- function(rhs)
#{
#    TOL <- 1e-07
#    n <- nrow(rhs)
#    Xj <- rhs
#    qrX <- qr(Xj, tol=TOL)
#    Q <- qr.Q(qrX)
#    H <- tcrossprod(Q[,1:qrX$rank])
#    return(H)
#}

#' INTERNAL: Calculates the hat matrix for the variable of interest or H2
#' 
#' @author Zarrar Shehzad
#' @param rhs right-hand matrix
#' @param grps vector indicating the factor index for each regressor
#' @param f.ind factor index to examine
#' @param o.inds vector of indices for observations (e.g., for a permutation)
#' @return matrix
mdmr_model.hat_matrix_2 <- function(rhs, grps, f.ind, o.inds=NULL)
{
    TOL <- 1e-07
    u.grps <- unique(grps)
    if (is.null(o.inds))
        o.inds <- 1:nrow(rhs)
    
    # H
    Xj <- rhs
    cols <- grps %in% u.grps[f.ind]
    Xj[,cols] <- Xj[o.inds,cols]
    qrX <- qr(Xj, tol=TOL)
    Q <- qr.Q(qrX)
    H <- tcrossprod(Q[,1:qrX$rank])
    
    # H2
    cols <- grps %in% u.grps[-f.ind]
    Xj <- rhs[,cols]
    qrX <- qr(Xj, tol = TOL)
    Q <- qr.Q(qrX)
    H2 <- H - tcrossprod(Q[, 1:qrX$rank])
    
    return(H2)
}

##' INTERNAL: Calculates the hat matrices for everything but the variable of 
##' interest or the H1s
##' 
##' @author Zarrar Shehzad
##' @param rhs Right-hand matrix
##' @param grps Vector indicating the factor index for each regressor
##' @param H Full-model hat matrix
##' @return list of matrices, one for each set of factors
#mdmr_model.hat_matrix_1 <- function(rhs, grps, H)
#{
#    TOL <- 1e-07
#    u.grps <- unique(grps)
#    
#    H1s <- lapply(1:length(u.grps), function(j) {
#        Xj <- rhs[, grps %in% u.grps[j]]
#        qrX <- qr(Xj, tol = TOL)
#        Q <- qr.Q(qrX)
#        H - tcrossprod(Q[, 1:qrX$rank])
#    })
#    
#    return(H1s)
#}

#' INTERNAL: Calculates 1 minus the full-model matrix (used for residual stuff)
#' 
#' Note: if o.inds is supplied, then the rows for the columns related to f.ind
#'       are permuted before calculated H and then I - H
#' 
#' @author Zarrar Shehzad
#' @param rhs right-hand matrix
#' @param grps vector indicating the factor index for each regressor
#' @param f.ind factor index to examine
#' @param o.inds vector of indices for observations (e.g., for a permutation)
#' @return matrix
mdmr_model.hat_matrix_ih <- function(rhs, grps, f.ind, o.inds=NULL)
{
    TOL <- 1e-07
    u.grps <- unique(grps)
    if (is.null(o.inds)) 
        o.inds <- 1:nrow(rhs)
    
    # H
    Xj <- rhs
    cols <- grps %in% u.grps[f.ind]
    Xj[,cols] <- Xj[o.inds,cols]
    qrX <- qr(Xj, tol=TOL)
    Q <- qr.Q(qrX)
    H <- tcrossprod(Q[,1:qrX$rank])
    
    # I - H
    IH <- diag(nrow(rhs)) - H
    
    return(IH)
}

#' Calculates the various hat matrices
#' 
#' @author Zarrar Shehzad
#' @param rhs Right-hand matrix
#' @param qrhs QR deconmposition of right-hand matrix
#' @return list of different hat matrices
mdmr_model.hat_matrices <- function(rhs, qrhs, perms=NULL)
{
    grps <- attr(qrhs, "grps")
    u.grps <- unique(grps)
    factors2perm <- attr(qrhs, "factors2perm")
    
    nobs <- nrow(rhs)
    if (is.null(perms))
        perms <- 1:nobs
    
    H2s <- lapply(factors2perm, function(j) 
                    mdmr_model.hat_matrix_2(rhs, grps, j, perms))
    
    IHs <- lapply(factors2perm, function(j) 
                    mdmr_model.hat_matrix_ih(rhs, grps, j, perms))
    
    # Return
    return(list(
        H2s = H2s, 
        IHs = IHs
    ))
}

#' Calculates the degrees of freedom
#' 
#' @author Zarrar Shehzad
#' @param qrhs QR decomposition of model matrix
#' @return list of dfs of residual and predictor variables
mdmr_model.df <- function(qrhs)
{
    n <- nrow(qrhs$qr)
    grps <- attr(qrhs, "grps")
    u.grps <- attr(qrhs, "u.grps")
    factors2perm <- attr(qrhs, "factors2perm")
    
    df.Res <- n - qrhs$rank
    df.Exp <- sapply(u.grps[factors2perm], function(i) sum(grps == i))
    
    return(list(
        res = df.Res, 
        exp = df.Exp
    ))
}



#parallel.perm = function(D.array, X, test.columns=ncol(X), permat=NULL, nperm=if (!is.null(permat)) ncol(permat) else 999, H2mat=NULL, IHmat=NULL, report.every=50) {
#
#	require(MASS)
#
#	n = nrow(X)
#
#	if (is.null(permat)) {
#
#		permat = matrix(NA, n, nperm)
#
#		for (k in 1:nperm) permat[ , k] = sample(n)
#
#	}
#
#	ndim = length(dim(D.array))
#
#    if (ndim==2) D.mat = D.array
#
#    else if (ndim==3) D.mat = matrix(D.array, ncol=dim(D.array)[3])
#
#    else stop("'D.array' must be either 3-way array of distance matrices, or a matrix whose columns are vectorized distance matrices")
#
#    if (nrow(D.mat) != n^2) stop("Distance matrices in 'D.array' must be n x n, where n is the number of rows of 'X'")
#    Amat = -D.mat[,]*D.mat[,]/2
#    H = X %*% ginv(X)
#    IH = diag(n) - H
#    H2 = H - X[ , -test.columns] %*% ginv(X[ , -test.columns])
#    if (is.null(H2mat)) {        
#
#        H2mat = IHmat = matrix(NA, n^2, nperm)
#        for (k in 1:nperm) {
#
#    	    if (k %% report.every==0) cat("Permutation", k, "\n")
#
#    	    prm = permat[ , k]
#
#    	    H2mat[ , k] = H2[prm, prm]
#
#    	    IHmat[ , k] = IH[prm, prm]
#
#    	}
#
#        # H2mat = H2mat - 1/n  # no!!!
#
#    }
#    F = ((n-ncol(X)) / length(test.columns)) * as.vector(crossprod(Amat, as.vector(H2)) / crossprod(Amat, as.vector(IH)))
#    permF = ((n-ncol(X)) / length(test.columns)) * crossprod(H2mat, Amat) / crossprod(IHmat, Amat)
#    maxpermF=apply(permF, 1, max)
#
#    list(F=F, permF=permF, pvalue=apply(rbind(F, permF), 2, function(v) sum(v[1]<=v)) / (nperm+1), maxpermF=maxpermF, corrected.p=(rowSums(outer(F, maxpermF, "<="))+1) / (nperm+1), H2mat=H2mat, IHmat=IHmat)
#
#}
#
#
#
#
#H2mat[ , k] = (I-H1) %*% H2[prm, prm] %*% (I-H1)
#IHmat[ , k] = (I-H1) %*% IH[prm, prm] %*% (I-H1)
