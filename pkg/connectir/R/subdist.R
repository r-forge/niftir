# This checks that everything in folder is good
check_subdist <- function(sdir) {
    
}

# This creates a new subdist directory and relevant files
# and returns a new subdist object
create_subdist <- function(outdir, infiles, masks, opts) {
    if (file.exists(outdir))
        stop("Output cannot exist")
    
    infuncdir <- file.path(outdir, "input_funcs")
    inmaskdir <- file.path(outdir, "input_masks")
    
    # Create directories
    dir.create(outdir)
    dir.create(infuncdir)
    dir.create(inmaskdir)
    
    # Create symlinks for the input funcs
    for (i in 1:length(infiles)) {
        from <- infiles[i]
        to <- file.path(infuncdir, sprintf("scan%04i.%s", i, getext(from)))
        file.symlink(from, to)
    }
    
    # Get a header file from the first functional
    header <- read.nifti.header(infiles[i])
    header$dim <- header$dim[1:3]
    header$pixdim <- header$pixdim[1:3]
    
    # Write the brain masks
    tmpnames <- names(masks)
    for (i in 1:length(masks)) {
        outfile <- file.path(inmaskdir, sprintf("%s.nii.gz", tmpnames[[i]]))
        write.nifti(masks[[i]], header, outfile=outfile)
    }
    
    # Copy over standard brains
    std2 <- system.file("data/standard_2mm.nii.gz", package="connectir")
    std4 <- system.file("data/standard_4mm.nii.gz", package="connectir")
    file.copy(std2, file.path(outdir, "standard_2mm.nii.gz"))
    file.copy(std4, file.path(outdir, "standard_4mm.nii.gz"))
    
    # Save options
    opts$outdir <- outdir
    opts$infiles <- infiles
    save(opts, file=file.path(outdir, "options.rda"))
    
    # Want to create subdist matrix
    nsubs <- length(infiles)
    nvoxs <- sum(masks$brainmask)
    big.matrix(nsubs^2, nvoxs, type="double")
}

compute_subdist <- function(funclist, subdist, seed_inds, blocksize, ztransform, start=1, verbose=TRUE) {
    nseeds <- length(seed_inds)
    blocks <- niftir.split.indices(start, nseeds, by=blocksize)
    
    dfun <- function(i, blocks, seed_inds, funclist, subdist, ztransform, verbose, pb) {
        if (verbose)
            update(pb, i)
        inds_CHUNK <- seed_inds[blocks$starts[i]:blocks$ends[i]]
        cormaps_list <- vbca_batch(funclist, inds_CHUNK, ztransform=ztransform)
        subdist_CHUNK <- sub.big.matrix(subdist, firstCol=blocks$starts[i], lastCol=blocks$ends[i])
        compute_subdist_worker(cormaps_list, inds_CHUNK, subdist_CHUNK)
    }
    
    if (verbose)
        pb <- progressbar(blocks$n)
    else
        pb <- NULL
    
    if (getDoParRegistered()) {
        foreach(i=1:blocks$n, .packages=c("connectir")) %dopar% 
            dfun(i, blocks, seed_inds, funclist, subdist, ztransform, verbose, pb)
    }
    else {
        foreach(i=1:blocks$n, .packages=c("connectir")) %do% 
            dfun(i, blocks, seed_inds, funclist, subdist, ztransform, verbose, pb)
    }
    
    if (verbose)
        end(pb)
}

compute_subdist_worker <- function(sub.cormaps, inds, outmat=NULL, type="double", ...) {
    nsubs <- length(sub.cormaps)
    nvoxs <- ncol(sub.cormaps[[1]])
    nseeds <- nrow(sub.cormaps[[1]])
    if (nseeds != length(inds))
        stop("length of inds doesn't match nrow of first sub.cormaps element")
    
    if (is.null(outmat))
        outmat <- big.matrix(nsubs^2, nseeds, type=type, ...)
    else if (ncol(outmat) != nseeds || nrow(outmat) != nsubs^2)
        stop("dimensions of outmat do not match nsubs and nseeds values")
    
    subsMap <- big.matrix(nvoxs-1, nsubs, type=type, ...)
    ALPHA <- 1/(nvoxs-2)
    voxs <- 1:nvoxs
    for (i in 1:nseeds) {
        .Call("CombineSubMapsMain", sub.cormaps, subsMap@address, as.double(i), as.double(voxs[-inds[i]]), as.double(nvoxs-1), as.double(nsubs))
        col <- sub.big.matrix(outmat, firstCol=i, lastCol=i)
        dgemm(C=col, A=subsMap, B=subsMap, TRANSA='t', ALPHA=ALPHA, LDC=as.double(nsubs))
        .Call("BigSubtractScalarMain", col@address, as.double(1), TRUE);
    }
    
    gc(F)
    
    return(outmat)
}

# bigmat: rows=subject distances, cols=voxels
# not that each column is a vectorized version of n x n matrix comparing subjects where n^2 = # of row elements
gower.subdist <- function(bigmat, gower.bigmat=NULL, verbose=TRUE, do.parallel=FALSE, ...) {
    # this all does
    # G <- -0.5 * -(dmat*dmat) %*% (I - ones %*% t(ones)/n)
    
    # setup
    nc <- ncol(bigmat)
    n <- sqrt(nrow(bigmat))
    I <- diag(n)
    ones <- matrix(1, nrow=n)
    if (is.null(gower.bigmat))
        gower.bigmat <- deepcopy(bigmat, ...)
    
    # bigmat <- bigmat * bigmat
    .Call("BigPowMain", gower.bigmat@address, as.double(2))
    
    # bigmat <- -(bigmat)/2
    dscal(ALPHA=-0.5, Y=gower.bigmat)
    
    # I - ones %*% t(ones)/n
    adj <- I - ones %*% t(ones)/n
    
    # newmat <- bigmat %*% adj
    if (verbose)
        pb <- progressbar(nc)
    for (i in 1:nc) {
        if (verbose)
            update(pb, i)
        gower.vox <- matrix(gower.bigmat[,i], n, n)
        gower.bigmat[,i] <- as.vector(gower.vox %*% adj)
    }
    if (verbose)
        end(pb)
    
    return(gower.bigmat)
}

gower.subdist <- function(bigmat, gower.bigmat=NULL, verbose=TRUE, ...) {
    # this all does
    # G <- -0.5 * -(dmat*dmat) %*% (I - ones %*% t(ones)/n)
    
    # setup
    nc <- ncol(bigmat)
    n <- sqrt(nrow(bigmat))
    I <- diag(n)
    ones <- matrix(1, nrow=n)
    if (is.null(gower.bigmat))
        gower.bigmat <- deepcopy(bigmat, ...)
    
    # bigmat <- bigmat * bigmat
    .Call("BigPowMain", gower.bigmat@address, as.double(2))
    
    # bigmat <- -(bigmat)/2
    dscal(ALPHA=-0.5, Y=gower.bigmat)
    
    # I - ones %*% t(ones)/n
    adj <- I - ones %*% t(ones)/n
    
    if (verbose)
        pb <- progressbar(nc)
    
    # newmat <- bigmat %*% adj
    gfun <- function(i, gower.bigmat) {
        if (verbose)
            update(pb, i)
        gower.vox <- matrix(gower.bigmat[,i], n, n)
        gower.bigmat[,i] <- as.vector(gower.vox %*% adj)
        return(NULL)
    }
    for (i in seq_len(nc))
        gfun(i, gower.bigmat)
    #TODO: seems like the parallel thing here takes longer than it should
    #if (getDoParWorkers() == 1) {
    #    for (i in seq_len(nc)) gfun(i, gower.bigmat)
    #} else {
    #    foreach(i = seq_len(nc)) %dopar% gfun(i, gower.bigmat)
    #}
    
    if (verbose)
        end(pb)
    
    return(gower.bigmat)
}

square.subdist <- function(bigmat, square.bigmat=NULL, ...) {
    # this all does
    # Amat <- -0.5 * -(dmat*dmat)
    
    # setup
    if (is.null(square.bigmat))
        square.bigmat <- deepcopy(bigmat, ...)
    
    # bigmat <- bigmat * bigmat
    .Call("BigPowMain", square.bigmat@address, as.double(2))
    
    # bigmat <- -(bigmat)/2
    dscal(ALPHA=-0.5, Y=square.bigmat)
    
    return(square.bigmat)
}
