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
        
        return(colmean(cormat_CHUNK))
    }
    
    ## initiate progress bar
    if (verbose)
        pb <- progressbar(blocks$n)
    else
        pb <- NULL
    
    ## loop
    if (getDoParRegistered() && getDoParWorkers() > 1) {
        lo <- min(getDoParWorkers()*3, blocks$n)
        superblocks <- niftir.split.indices(1, blocks$n, length.out=lo)
        gcor <- foreach(si=1:superblocks$n, .packages=c("connectir"), .combine='c') %dopar% 
            sapply(superblocks$starts[si]:superblocks$ends[si], dfun)
    }
    else {
        gcor <- sapply(1:blocks$n, dfun)
    }
    
    ## end progress bar
    if (verbose)
        end(pb)
    
    return(gcor)
}

# Only for testing purposes
.kendall <- function(ratings) {
    ratings <- as.matrix(na.omit(ratings))
    ns <- nrow(ratings)
    nr <- ncol(ratings)
    
    ratings.rank <- apply(ratings, 2, rank)
    coeff <- (12 * var(apply(ratings.rank, 1, sum)) * (ns - 
        1))/(nr^2 * (ns^3 - ns))
        
    return(coeff)
}

reho_worker <- function(bigmat, inds, ...) {
    return(bigmat[,inds])
}

reho <- function(bigmat, header, mask, nei=1, min.nei=0.25, verbose=TRUE, FUN=reho_worker, ...) {
    header$dim <- header$dim[1:3]
    header$pixdim <- header$pixdim[1:3]
    xcoords <- coords(header, mask)
    nvoxs <- ncol(bigmat)
    if (sum(mask) != nvoxs)
        stop("Number of TRUE elements in mask does not equal number of columns in bigmat")
    
    relative_coords <- expand.grid(list(i=-nei:nei, j=-nei:nei, k=-nei:nei))
    
    opts <- list()
    opts$min.nei <- ceiling(nrow(relative_coords)*min.nei)
    opts$min.i <- 0
    opts$max.i <- nvoxs+1
    opts$dim <- header$dim
    
    arr2vec <- function(df, dim) {
        (df$k-1)*dim[1]*dim[2] + (df$j-1)*dim[1] + df$i
    }
    
    rfun <- function(i, ...) {
        if (verbose)
            update(pb, i)
        
        inds <- arr2vec(xcoords[i,] + relative_coords, opts$dim)
        touse <- inds>opts$min.i & inds<opts$max.i
        if (sum(touse) < opts$min.nei)
            return(0)
        inds <- inds[touse]
        
        return(FUN(bigmat, inds, ...))
    }
    
    if (verbose)
        pb <- progressbar(nvoxs)
    else
        pb <- NULL
    
    if (getDoParRegistered() && getDoParWorkers() > 1) {
        blocks <- niftir.split.indices(1, nvoxs, length.out=getDoParWorkers())
        reho.vals <- foreach(bi=1:blocks$n, .inorder=TRUE, .combine='c') %dopar% 
            sapply(blocks$starts[bi]:blocks$ends[bi], rfun, ...)
    }
    else {
        reho.vals <- sapply(1:nvoxs, rfun, ...)
    }
    
    if (verbose)
        end(pb)
    
    return(reho.vals)   # todo: convert to nifti object?
}

