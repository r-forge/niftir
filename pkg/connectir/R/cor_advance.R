gcor <- function(bigmat, blocksize, ztransform=FALSE, verbose=TRUE) {
    ## input info
    nvoxs <- ncol(bigmat)
    vox_inds <- 1:nvoxs
    blocks <- niftir.split.indices(1, nvoxs, by=blocksize)
    
    ## function that computes correlation maps for subset of voxels
    ## + gets average of each correlation map
    dfun <- function(i) {
        if (verbose)
            update(pb, i)
        
        inds_CHUNK <- vox_inds[blocks$starts[i]:blocks$ends[i]]
        cormat_CHUNK <- vbca(bigmat, inds_CHUNK, ztransform=ztransform)
        
        return(bm_rowmean(cormat_CHUNK))
    }
    
    ## initiate progress bar
    if (verbose)
        pb <- progressbar(blocks$n)
    
    ## loop
    if (getDoParRegistered() && getDoParWorkers() > 1) {
        lo <- min(getDoParWorkers()*3, blocks$n)
        superblocks <- niftir.split.indices(1, blocks$n, length.out=lo)
        gcor <- foreach(si=1:superblocks$n, .packages=c("connectir")) %dopar% 
            unlist(lapply(superblocks$starts[si]:superblocks$ends[si], dfun))
    }
    else {
        gcor <- lapply(1:blocks$n, dfun)
    }
    
    gcor <- unlist(gcor)
    
    ## end progress bar
    if (verbose)
        end(pb)
    
    return(gcor)
}

# This is a simple version of kendall
kendall_ref <- function(ratings) {
    ratings <- as.matrix(na.omit(ratings))
    ns <- nrow(ratings)
    nr <- ncol(ratings)
    
    ratings.rank <- apply(ratings, 2, rank)
    coeff <- (12 * var(apply(ratings.rank, 1, sum)) * (ns - 
        1))/(nr^2 * (ns^3 - ns))
        
    return(coeff)
}

# This computes a kendall's W examining the consistency of
# each voxel's connectivity map across participants
kendall <- function(subs.bigmats, blocksize, ztransform=FALSE, verbose=TRUE) {
    ## input info
    nsubs <- length(subs.bigmats)
    nvoxs <- ncol(subs.bigmats[[1]])
    vox_inds <- 1:nvoxs
    blocks <- niftir.split.indices(1, nvoxs, by=blocksize)
    
    ## initiate progress bar
    if (verbose)
        pb <- progressbar(blocks$n)
    
    ## function
    kfun_worker <- function(ratings) {
        var(rowSums(apply(ratings, 2, rank)))
    }
    
    kfun <- function(i, ns, nr) {
        cols <- vox_inds[blocks$starts[i]:blocks$ends[i]]
        ncols <- length(cols)
        
        cormats <- vbca_batch(subs.bigmats, cols, ztransform)
        
        vals <- sapply(1:ncols, function(ci) {
            kfun_worker(sapply(cormats, function(mat) mat[ci,]))
        })
        coeffs <- (12 * vals * (ns-1))/(nr^2 * (ns^3 - ns))
        
        if (verbose)
            update(pb, i)
        
        return(coeffs)
    }
    
    if (getDoParRegistered() && getDoParWorkers() > 1) {
        lo <- min(getDoParWorkers()*3, blocks$n)
        superblocks <- niftir.split.indices(1, blocks$n, length.out=lo)
        gcor <- foreach(si=1:superblocks$n, .packages=c("connectir")) %dopar% 
            unlist(lapply(superblocks$starts[si]:superblocks$ends[si], kfun, nvoxs, nsubs))
        gcor <- unlist(gcor)
    }
    else {
        gcor <- unlist(lapply(1:blocks$n, kfun, nvoxs, nsubs))
    }
    
    if (verbose)
        end(pb)
    
    return(gcor)
}

reho_worker <- function(bigmat, inds, ...) {
    return(kendall(bigmat[,inds], ...)$value)
}

reho <- function(bigmat, nei=1, nei.dist=3, min.nei=0.5, verbose=TRUE, FUN=reho_worker, ...) {
    header <- bigmat@header
    header$dim <- header$dim[1:3]
    header$pixdim <- header$pixdim[1:3]
    
    mask <- bigmat@mask
    nvoxs <- ncol(bigmat)
    if (sum(mask) != nvoxs)
        stop("Number of TRUE elements in mask does not equal number of columns in bigmat")
    
    dims <- header$dim
    moffsets <- expand.grid(list(i=-nei:nei, j=-nei:nei, k=-nei:nei))
    dist <- rowSums(abs(moffsets))
    moffsets <- moffsets[dist<=nei.dist,]
    offsets <- moffsets$k*dims[1]*dims[2] + moffsets$j*dims[1] + moffsets$i
    
    mat2arr_inds <- which(mask)
    arr2mat_inds <- mask*1
    arr2mat_inds[mask] <- 1:nvoxs
    
    opts <- list()
    opts$min.nei <- ceiling(nrow(moffsets)*min.nei)
    opts$rawi.min <- 0
    opts$rawi.max <- length(mask)+1
    opts$dim <- header$dim
    opts$mat2arr_inds <- mat2arr_inds
    opts$arr2mat_inds <- arr2mat_inds
    opts$offsets <- offsets
    
    # have list of mask inds
    
    rfun <- function(i, bm, opts, ...) {
        if (verbose)
            update(pb, i)
        
        raw_inds <- opts$offsets + opts$mat2arr_inds[i]
        raw_inds <- raw_inds[raw_inds > opts$rawi.min & raw_inds < opts$rawi.max]
        filt_inds <- opts$arr2mat_inds[raw_inds]
        filt_inds <- filt_inds[filt_inds>0]
        
        if (length(filt_inds) < opts$min.nei)
            return(0)
        
        return(FUN(bm[,filt_inds], ...))        
    }
    
    if (verbose)
        pb <- progressbar(nvoxs)
    else
        pb <- NULL
    
    if (getDoParRegistered() && getDoParWorkers() > 1) {
        blocks <- niftir.split.indices(1, nvoxs, length.out=getDoParWorkers())
        reho.vals <- foreach(bi=1:blocks$n, .inorder=TRUE) %dopar% 
            sapply(blocks$starts[bi]:blocks$ends[bi], rfun, bigmat, opts, ...)
        reho.vals <- unlist(reho.vals)
    }
    else {
        reho.vals <- sapply(1:nvoxs, rfun, bigmat, opts, ...)
    }
    
    if (verbose)
        end(pb)
    
    return(list(reho=reho.vals, hdr=header, mask=mask))
}

