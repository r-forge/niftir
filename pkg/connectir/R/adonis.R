# BELOW: CODE DIRECTLY TAKEN FROM vegan
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
# ABOVE: CODE DIRECTLY TAKEN FROM vegan

# NOTE: assume that y-intercept exists...need to deal with case when don't have intercept
mdmr.prepare.model <- function(formula, model, contr.unordered="contr.sum", contr.ordered="contr.poly") {
    n <- nrow(model)
    Terms <- terms(formula, data = model)
    formula[[2]] <- NULL
    
    # Get right hand matrix
    rhs.frame <- model.frame(formula, model, drop.unused.levels = TRUE)
    op.c <- options()$contrasts
    options(contrasts = c(contr.unordered, contr.ordered))
    rhs <- model.matrix(formula, rhs.frame)
    options(contrasts = op.c)
    grps <- attr(rhs, "assign")
    qrhs <- qr(rhs)
    rhs <- rhs[, qrhs$pivot, drop = FALSE]
    rhs <- rhs[, 1:qrhs$rank, drop = FALSE]
    
    # Check if rank deficient
    solve(t(rhs) %*% rhs)
    
    # Get different groups of variables
    grps <- grps[qrhs$pivot][1:qrhs$rank]
    u.grps <- unique(grps)
    nterms <- length(u.grps) - 1
    factor.names <- attr(attr(rhs.frame, "terms"), "term.labels")[u.grps]
    
    # Get hat matrix
    H <- rhs %*% ginv(rhs)
    IH <- diag(n) - H
    
    # Create individual hat matrices that resolve each factor
    H2s <- lapply(2:length(u.grps), function(j) {
        iis <- grps!=u.grps[j]
        H - rhs[,iis] %*% ginv(rhs[,iis])
    })
    
    # Get degrees of freedom
    df.Res <- n - qrhs$rank
    df.Exp <- sapply(u.grps[-1], function(i) sum(grps == i))
    
    return(list(
        rhs = rhs,
        grps = grps,
        u.grps = u.grps,
        factor.names = factor.names,
        H = H,
        H2s = H2s,
        IH = IH,
        df.Res = df.Res,
        df.Exp = df.Exp,
        nobs = n,
        nfactors = length(factor.names)
    ))
}

mdmr.prepare.permutations <- function(modelinfo, nperms, strata, max.iter, factors2perm) {
    # get matrix containing different permutations of observations
    p <- sapply(1:nperms, function(x) permuted.index(modelinfo$nobs, strata))
    # remove permutations that are significantly correlated with your model
    if (max.iter > 0) {
        rthresh <- tanh(1.96/sqrt(modelinfo$nobs-3)) # r value that is significant (p < 0.05)
        which.cor <- 1:nperms
        for (iter in 1:max.iter) {
            which.cor <- unique(unlist(lapply(factors2perm, function(i) {
                H <- modelinfo$H2s[[i]]
                tmp <- sapply(which.cor, function(j) cor(H[,1], H[p[,j],1]))
                which(abs(tmp)>rthresh)
            })))
            if (length(which.cor)==0)
                break
            p[,which.cor] <- sapply(1:length(which.cor), function(x) permuted.index(modelinfo$nobs, strata))
        }
    }
    return(p)
}

mdmr.prepare.permH2s <- function(modelinfo, p) {
    nperms <- ncol(p) + 1   # original data + nperms
    nfactors <- length(modelinfo$factors2perm)
    H2mats <- lapply(1:nfactors, function(ii) {
        i <- modelinfo$factors2perm[ii]
        bigmat <- big.matrix(modelinfo$nobs^2, nperms, type="double")
        bigmat[,1] <- as.vector(modelinfo$H2s[[i]])
        for (k in 2:nperms) {
            prm <- p[,k-1]
            bigmat[,k] <- as.vector(modelinfo$H2s[[i]][prm,prm])
        }
        return(bigmat)
    })
    return(H2mats)
}

mdmr.prepare.permIH <- function(modelinfo, p) {
    nperms <- ncol(p) + 1   # original data + nperms
    IHmat <- big.matrix(modelinfo$nobs^2, nperms, type="double")
    IHmat[,1] <- as.vector(modelinfo$IH)
    for (k in 2:nperms) {
        prm <- p[,k-1]
        IHmat[,k] <- as.vector(modelinfo$IH[prm,prm])
    }
    return(IHmat)
}

mdmr_worker <- function(firstVox, lastVox, Gmat, H2mats, IHmat, df.Res, df.Exp, Pmat, Fperms) {
    # n = number of subjects
    # where H2mat has rows of H's (n^2 rows) and cols of # of permutations
    # where Gmat has rows of dmats (n^2 rows)  and cols of # of voxels
    # result has rows of # of permutations and cols of # of voxels
    
    inds <- firstVox:lastVox
    nperms <- ncol(IHmat)
    nvoxs <- length(inds)
    nterms <- length(H2mats)
    
    # get part of bigmatrix with Gower's centered subject distances
    Gmat.chunk <- deepcopy(Gmat, cols=inds)
    
    # get part of Pmat
    Pmat.chunk <- sub.big.matrix(Pmat, firstRow=firstVox, lastRow=lastVox)
    
    # compute error term
    error.variance <- big.matrix(nperms, nvoxs, type="double")
    dgemm(C=error.variance, A=IHmat, B=Gmat.chunk, TRANSA='t')
    
    explained.variance <- big.matrix(nperms, nvoxs, type="double")
    for (i in 1:nterms) {
        # compute explained variance
        dgemm(C=explained.variance, A=H2mats[[i]], B=Gmat.chunk, TRANSA='t')
        
        # explained-variance / error
        Fstats <- sub.big.matrix(Fperms[[i]], firstCol=firstVox, lastCol=lastVox)
        do.operator(explained.variance, error.variance, "/", Fstats)
        
        # adjust for degrees of freedom
        dscal(Y=Fstats, ALPHA=df.Res/df.Exp[i])
        
        # get pvals (TODO: convert below line to C++ code)
        #Pmat.chunk[,i] <- apply(Fstats[,], 2, function(x) sum(x >= x[1])/nperms)
        .Call("ComputePvalsMain", Fstats@address, Pmat.chunk@address, as.double(i), PACKAGE="connectir")
    }
    
    # cleanup
    rm(explained.variance)
    rm(error.variance)
    gc(FALSE)
    
    return(NULL)
}


# assume each column of x has been gower centered
mdmr <- function(x, formula, model, nperms=4999, factors2perm=NULL, voxs=1:ncol(x), block.size=250, verbose=TRUE, contr.unordered="contr.sum", contr.ordered="contr.poly", max.iter=10, strata=NULL) {
    # todo: test if x is matrix or big.matrix?
    
    if (!is.data.frame(model))
        stop("'model' input must be a data frame")
    
    nVoxs <- length(voxs)
    nSubs <- sqrt(nrow(x))  # don't want length(subs) here
    blocks <- niftir.split.indices(1, nVoxs, by=block.size)
    
    # Prepare model matrix
    vcat(verbose, "Preparing model stuff")
    modelinfo <- mdmr.prepare.model(formula, model, contr.unordered, contr.ordered)
    
    # Permutation Business
    vcat(verbose, "Preparing permutation related inputs")
    ## factors to permute
    if (is.null(factors2perm))
        factors2perm <- 1:modelinfo$nfactors
    else if (is.character(factors2perm))
        factors2perm <- sapply(factors2perm, function(x) which(modelinfo$factor.names==x))
    factor.names <- modelinfo$factor.names[factors2perm]
    ## get permutations
    p <- mdmr.prepare.permutations(modelinfo, nperms, strata, max.iter, factors2perm)
    ## refine df of experiment
    modelinfo$df.Exp <- sapply(factors2perm, function(i) modelinfo$df.Exp[i])
    modelinfo$factors2perm <- factors2perm
    
    # Hat matrices prepared with all possible permuted models
    vcat(verbose, "Preparing permuted model matrices")
    H2mats <- mdmr.prepare.permH2s(modelinfo, p)
    IHmat <- mdmr.prepare.permIH(modelinfo, p)
    
    # Create output matrices
    vcat(verbose, "Preparing output matrices (pvals and fstats)")
    ## pvals
    Pmat <- big.matrix(nVoxs, length(factors2perm), type="double")
    ## fstats
    Fperms <- lapply(1:length(factors2perm), function(i) {
        big.matrix(nperms+1, nVoxs, type="double")
    })
    
    vcat(verbose, "Will calculate permutation based p-values for the following factors:", factor.names)
    
    vcat(verbose, "Computing MDMR across", blocks$n, "blocks")
    ## progress bar
    if (verbose)
        pb <- progressbar(blocks$n)
    ## loop through
    if (getDoParRegistered() && getDoParWorkers() > 1) {
        foreach(i=1:blocks$n) %dopar% {
            if (verbose)
                update(pb, i)
            mdmr_worker(blocks$starts[i], blocks$ends[i], x, H2mats, IHmat, modelinfo$df.Res, modelinfo$df.Exp, Pmat, Fperms)
        }
    } else {
        for (i in 1:blocks$n) {
            if (verbose)
                update(pb, i)
            mdmr_worker(blocks$starts[i], blocks$ends[i], x, H2mats, IHmat, modelinfo$df.Res, modelinfo$df.Exp, Pmat, Fperms)
        }
    }
    ## end progress bar
    if (verbose)
        end(pb)
    
    structure(
        list(
            modelinfo=modelinfo,
            pvals=Pmat,
            fstats=Fperms,
            perms=p
        ),
        class="mdmr"
    )
}

save_mdmr <- function(obj, sdir, mdir, formula, verbose=TRUE) {
    if (!file.exists(sdir))
        stop("Cannot save MDMR to ", sdir, " since it doesn't exist")
    
    vcat(verbose, "...creating MDMR output directory")
    mdmr.output <- mdir
    if (file.exists(mdmr.output))
        stop("MDMR output ", mdmr.output, " already exists")
    dir.create(mdmr.output)
    
    mpath <- function(...) file.path(mdmr.output, ...)
    
    vcat(verbose, "...reading in brain mask")
    seedfn <- file.path(sdir, "input_masks", "seedmask.nii.gz")
    mask <- read.mask(seedfn)
    header <- read.nifti.header(seedfn)
    
    # Model Stuff
    modelinfo <- obj$modelinfo
    modelinfo$perms <- obj$perms
    ## save
    vcat(verbose, "...saving model info")
    save(modelinfo, file=mpath("modelinfo.rda"))
    ## formula
    cat(as.character(formula)[-c(1:2)], file=mpath("formula.txt"))
    ## factors
    vcat(verbose, "...saving factor names")
    factornames <- modelinfo$factor.names[modelinfo$factors2perm]
    nfactors <- length(factornames)
    totext <- sprintf("# All Factors\n%s\n# Permuted Factors\n%s\n", 
                modelinfo$factor.names, factornames)
    cat(totext, file=mpath("factorinfo.txt"))
    ## evs
    vcat(verbose, "...saving evs")
    evs <- modelinfo$rhs
    write.csv(evs, file=mpath("evs.csv"))
    ## models
    # H2s <- modelinfo$H2s
    ## residuals
    vcat(verbose, "...saving residuals")
    IH <- modelinfo$IH
    write.table(IH[,], quote=F, row.names=F, col.names=F, file=mpath("residuals.2D"))
    
    vcat(verbose, "...saving p-values")
    Pmat <- obj$pvals
    for (i in 1:nfactors) {
        fn <- mpath(sprintf("pvals_%s.nii.gz", factornames[i]))
        write.nifti(Pmat[,i], header, mask, outfile=fn)
    }
    
    vcat(verbose, "...saving FDR corrected p-values")
    CorrPmat <- big.matrix(nrow(Pmat), ncol(Pmat), type="double")
    for (i in 1:nfactors) {
        fn <- mpath(sprintf("fdr_pvals_%s.nii.gz", factornames[i]))
        tmp <- p.adjust(Pmat[,i], "BH")
        CorrPmat[,i] <- tmp
        write.nifti(tmp, header, mask, outfile=fn)
        rm(tmp)
        gc(FALSE)
    }
    
    vcat(verbose, "...saving z-statics")
    for (i in 1:nfactors) {
        fn <- mpath(sprintf("zstats_%s.nii.gz", factornames[i]))
        tmp <- qt(Pmat[,i], Inf, lower.tail=FALSE)
        write.nifti(tmp, header, mask, outfile=fn)
        rm(tmp)
        gc(FALSE)
    }
    
    vcat(verbose, "...saving FDR corrected z-statistics")
    for (i in 1:nfactors) {
        fn <- mpath(sprintf("fdr_zstats_%s.nii.gz", factornames[i]))
        tmp <- qt(CorrPmat[,i], Inf, lower.tail=FALSE)
        write.nifti(tmp, header, mask, outfile=fn)
        rm(tmp)
        gc(FALSE)
    }
    rm(CorrPmat)
    gc(FALSE)
    
    vcat(verbose, "...saving permutated F-statistics")
    Fperms <- obj$fstats
    for (i in 1:nfactors) {
        tmp <- deepcopy(Fperms[[i]], backingpath=mdmr.output, 
            backingfile=sprintf("fperms_%s.bin", factornames[i]), 
            descriptorfile=sprintf("fperms_%s.desc", factornames[i]))
        rm(tmp)
        gc(FALSE)
    }
}


