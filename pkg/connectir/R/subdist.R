check_dmat <- function(dmat) {
    TOL <- .Machine$double.eps ^ 0.5
    dmat <- abs(as.matrix(dmat))
    diag_dmat <- diag(dmat)
    off_dmat <- dmat[lower.tri(dmat)]
    
    if (any(is.na(dmat)))
        stop("NAs were present in distance matrix")
    
    if (all(dmat < TOL))
        stop("All zeros in distance matrix")
    
    if (all(diag_dmat>TOL)) {
        warning("Diagonal of distance matrix is all non-zeros\n")
    } else if (any(diag_dmat>TOL)) {
        cat("Diagonal of distance matrix has non-zeros\n")
        print(diag_dmat)
    }
    
    if (any(off_dmat<TOL))
        cat("Off-diagonal of distance matrix has some zeros\n")
}

check_gmat <- function(gmat) {
    TOL <- .Machine$double.eps ^ 0.5
    dmat <- abs(as.matrix(gmat))
    
    if (any(is.na(gmat)))
        stop("NAs were present in gower's centered distance matrix")
    if (all(dmat < TOL))
        stop("All zeros in gower's centered distance matrix")
}

# This checks that everything in folder is good
check_subdist <- function(sdir) {
    # Checking Functions
    checkpath <- function(p, is.dir) { 
        if (!file.exists(p)) {
            if (is.dir)
                stop("Directory: ", p, " does not exist but required")
            else
                stop("File: ", p, " does not exist but required")
        }
    }
    checkthing <- function(comparison, ...) {
        if (!comparison)
            stop(...)
    }
    
    # File/directory paths to check
    sdir <- abspath(sdir)
    optsfile <- file.path(sdir, "options.rda")
    infuncdir <- file.path(sdir, "input_funcs")
    inmaskdir <- file.path(sdir, "input_masks")
    inmaskfiles <- file.path(inmaskdir, 
        c("brainmask.nii.gz",  "seedmask.nii.gz")
    )
    seedfile <- file.path(inmaskdir, "seedmask.nii.gz")
    distfiles <- file.path(sdir, 
        c("subdist.desc", "subdist.bin", "subdist_gower.desc", "subdist_gower.bin")
    )
    sdistfile <- file.path(sdir, "subdist.desc")
    gdistfile <- file.path(sdir, "subdist_gower.desc")
    
    # Check main directory
    checkpath(sdir, TRUE)
    
    # Check opts and read them in
    checkpath(optsfile, FALSE)
    opts <- NULL
    load(optsfile)
    if (is.null(opts))
        stop("Optsfile doesn't have opts variable!")
    
    # Check input funcs
    checkpath(infuncdir, TRUE)
    checkthing(
        length(list.files(infuncdir))==length(opts$infiles), 
        "Missing some input functional files"
    )
    
    # Check input masks
    checkpath(inmaskdir, TRUE)
    lapply(inmaskfiles, checkpath, FALSE)
    
    # Check subdist and related files
    lapply(distfiles, checkpath, FALSE)
    
    # Count number of seed voxels
    mask <- read.mask(seedfile)
    nvoxs <- sum(mask)
    
    # Read in subdist and seedmask + check
    lapply(c(sdistfile, gdistfile), function(f) {
        x <- attach.big.matrix(f)
        ## make sure appropriate number of voxels
        checkthing(
            nvoxs==ncol(x),
            "# of seed voxels does not match # of columns in ", f
        )
        ## mask sure appropriate number of subjects
        nsubs <- sqrt(nrow(x))
        checkthing(
            length(opts$infiles)==nsubs,
            "# of subjects does not match # of rows in ", f
        )
    })
}

# This creates a new subdist directory and relevant files
# and returns a new subdist object
create_subdist <- function(outdir, infiles, masks, opts, ...) {
    if (file.exists(outdir))
        stop("Output cannot exist")
    
    infuncdir <- file.path(outdir, "input_funcs")
    inmaskdir <- file.path(outdir, "input_masks")
    
    # Create directories
    dir.create(outdir)
    dir.create(infuncdir)
    dir.create(inmaskdir)
    
    # Create symlinks for the input funcs
    if (!opts$"no-link-functionals") {
        for (i in 1:length(infiles)) {
            from <- infiles[i]
            to <- file.path(infuncdir, sprintf("scan%04i.%s", i, getext(from)))
            file.symlink(from, to)
        }
    }
    
    # Get a header file from the first functional
    header <- read.nifti.header(infiles[1])
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
    big.matrix(nsubs^2, nvoxs, type="double", ...)
}

compute_subdist <- function(funclist, subdist, seed_inds, blocksize, ztransform, start=1, verbose=TRUE, testonly=FALSE) {
    nseeds <- length(seed_inds)
    blocks <- niftir.split.indices(start, nseeds, by=blocksize)
    
#    dfun <- function(i, blocks, seed_inds, funclist, subdist, ztransform, verbose, pb) {
    dfun <- function(i, ...) {
        if (verbose) {
            update(pb, i)
            #msg <- sprintf("\nblock %i with voxels %i:%i\n", i, blocks$starts[i], blocks$ends[i])
            #cat(msg)
        }
        dist_inds_CHUNK <- blocks$starts[i]:blocks$ends[i]
        scor_inds_CHUNK <- seed_inds[dist_inds_CHUNK]
        cormaps_list <- vbca_batch(funclist, scor_inds_CHUNK, ztransform=ztransform, shared=FALSE)        
        tmp <- compute_subdist_worker(cormaps_list, scor_inds_CHUNK, subdist, dist_inds_CHUNK)
        rm(dist_inds_CHUNK, scor_inds_CHUNK, cormaps_list, tmp)
        gc(FALSE, TRUE)
        return(NULL)
    }
    
    # Test
    i <- 1
    if (verbose) {
        cat("...running a test (", blocks$starts[i],  ")\n")
        pb <- progressbar(i)
    } else {
        pb <- NULL
    }
    dfun(i)
    check_dmat(matrix(subdist[,blocks$starts[i]], sqrt(nrow(subdist))))
    check_dmat(matrix(subdist[,blocks$ends[i]], sqrt(nrow(subdist))))
    if (verbose)
        end(pb)
    if (testonly) {
        cat("...test only...\n")
        return(NULL)
    }
    
    # Subdist Calculation
    if (verbose) {
        cat("...now the real deal\n")
        pb <- progressbar(blocks$n)
    } else {
        pb <- NULL
    }
    
    if (getDoParRegistered() && getDoParWorkers() > 1) {
        lo <- min(getDoParWorkers()*3, blocks$n-1)
        superblocks <- niftir.split.indices(2, blocks$n, length.out=lo)
        foreach(si=1:superblocks$n, .packages=c("connectir"), .inorder=TRUE) %dopar% 
            for(i in superblocks$starts[si]:superblocks$ends[si]) dfun(i)
    }
    else {
        for (i in 2:blocks$n)
            dfun(i)
    }
    
    if (verbose)
        end(pb)
}

compute_subdist_worker <- function(sub.cormaps, cor_inds, outmat, dist_inds, type="double", ...) {
    nsubs <- length(sub.cormaps)
    nvoxs <- ncol(sub.cormaps[[1]])
    nseeds <- nrow(sub.cormaps[[1]])
    if (nseeds != length(cor_inds) || nseeds != length(dist_inds))
        stop("length of inds doesn't match nrow of first sub.cormaps element")
    
    #if (is.null(outmat))
    #    outmat <- big.matrix(nsubs^2, nseeds, type=type, ...)
    #else if (ncol(outmat) != nseeds || nrow(outmat) != nsubs^2)
    #    stop("dimensions of outmat do not match nsubs and nseeds values")
    
    subsMap <- big.matrix(nvoxs-1, nsubs, type=type, shared=FALSE, ...)
    ALPHA <- 1/(nvoxs-2)
    voxs <- 1:nvoxs
    for (i in 1:nseeds) {
        .Call("CombineSubMapsMain", sub.cormaps, subsMap@address, as.double(i), as.double(voxs[-cor_inds[i]]), as.double(nvoxs-1), as.double(nsubs))
        col <- sub.big.matrix(outmat, firstCol=dist_inds[i], lastCol=dist_inds[i])
        dgemm(C=col, A=subsMap, B=subsMap, TRANSA='t', ALPHA=ALPHA, LDC=as.double(nsubs))
        .Call("BigSubtractScalarMain", col@address, as.double(1), TRUE);
    }
    
    rm(subsMap)
    gc(F)
    
    return(outmat)
}

# bigmat: rows=subject distances, cols=voxels
# not that each column is a vectorized version of n x n matrix comparing subjects where n^2 = # of row elements
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
    gfun <- function(i) {
        if (verbose)
            update(pb, i)
        gower.vox <- matrix(gower.bigmat[,i], n, n)
        gower.bigmat[,i] <- as.vector(gower.vox %*% adj)
        return(NULL)
    }
    for (i in seq_len(nc))
        gfun(i)
    #TODO: seems like the parallel thing here takes longer than it should
    #if (getDoParWorkers() == 1) {
    #    for (i in seq_len(nc)) gfun(i)
    #} else {
    #    foreach(i = seq_len(nc)) %dopar% gfun(i)
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

slice.subdist <- function(bigmat, subs=1:sqrt(nrow(bigmat)), voxs=1:ncol(bigmat), ...) {
    matinds <- matrix(1:nrow(bigmat), sqrt(nrow(bigmat)))
    matinds <- as.vector(matinds[subs,subs])
    deepcopy(bigmat, cols=voxs, rows=matinds, ...)
}




### NEW CODE USING ARMADILLO

compute_subdist2 <- function(funclist, subdist, seed_inds, blocksize, ztransform, start=1, verbose=TRUE, testonly=FALSE, design_matrix=NULL) {
    nseeds <- length(seed_inds)
    blocks <- niftir.split.indices(start, nseeds, by=blocksize)
    
#    dfun <- function(i, blocks, seed_inds, funclist, subdist, ztransform, verbose, pb) {
    if (is.null(design_matrix)) {
        dfun <- function(i, ...) {
            if (verbose) {
                update(pb, i)
                #msg <- sprintf("\nblock %i with voxels %i:%i\n", i, blocks$starts[i], blocks$ends[i])
                #cat(msg)
            }
            
            dist_inds_CHUNK <- blocks$starts[i]:blocks$ends[i]
            scor_inds_CHUNK <- seed_inds[dist_inds_CHUNK]
            cormaps_list <- vbca_batch2(funclist, scor_inds_CHUNK, ztransform=ztransform, shared=FALSE)        
            tmp <- compute_subdist_worker2(cormaps_list, scor_inds_CHUNK, subdist, dist_inds_CHUNK)
            
            rm(dist_inds_CHUNK, scor_inds_CHUNK, cormaps_list, tmp)
            gc(FALSE)
            return(NULL)
        }
    } else {
        dfun <- function(i, ...) {
            if (verbose) {
                update(pb, i)
                #msg <- sprintf("\nblock %i with voxels %i:%i\n", i, blocks$starts[i], blocks$ends[i])
                #cat(msg)
            }
            dist_inds_CHUNK <- blocks$starts[i]:blocks$ends[i]
            scor_inds_CHUNK <- seed_inds[dist_inds_CHUNK]
            cormaps_list <- vbca_batch2(funclist, scor_inds_CHUNK, ztransform=ztransform, 
                                        shared=FALSE)        
            tmp <- compute_subdist_worker2_regress(cormaps_list, scor_inds_CHUNK, subdist, 
                                                   dist_inds_CHUNK, design_matrix)
            
            rm(dist_inds_CHUNK, scor_inds_CHUNK, cormaps_list, tmp)
            gc(FALSE)
            return(NULL)
        }
    }
    
    # Test
    i <- 1
    if (verbose) {
        cat("...running a test (", blocks$starts[i],  ")\n")
        pb <- progressbar(i)
    } else {
        pb <- NULL
    }
    dfun(i)
    check_dmat(matrix(subdist[,blocks$starts[i]], sqrt(nrow(subdist))))
    check_dmat(matrix(subdist[,blocks$ends[i]], sqrt(nrow(subdist))))
    if (verbose)
        end(pb)
    if (testonly) {
        cat("...test only...\n")
        return(NULL)
    }
    
    # Subdist Calculation
    if (verbose) {
        cat("...now the real deal\n")
        pb <- progressbar(blocks$n)
    } else {
        pb <- NULL
    }
    
    if (getDoParRegistered() && getDoParWorkers() > 1) {
        lo <- min(getDoParWorkers()*3, blocks$n-1)
        superblocks <- niftir.split.indices(2, blocks$n, length.out=lo)
        foreach(si=1:superblocks$n, .packages=c("connectir"), .inorder=TRUE) %dopar% 
            for(i in superblocks$starts[si]:superblocks$ends[si]) dfun(i)
    }
    else {
        for (i in 2:blocks$n)
            dfun(i)
    }
    
    if (verbose)
        end(pb)
}

compute_subdist_worker2 <- function(sub.cormaps, cor_inds, outmat, dist_inds, type="double", ...) {
    nsubs <- length(sub.cormaps)
    nvoxs <- ncol(sub.cormaps[[1]])
    nseeds <- nrow(sub.cormaps[[1]])
    if (nseeds != length(cor_inds) || nseeds != length(dist_inds))
        stop("length of inds doesn't match nrow of first sub.cormaps element")
    
    #if (is.null(outmat))
    #    outmat <- big.matrix(nsubs^2, nseeds, type=type, ...)
    #else if (ncol(outmat) != nseeds || nrow(outmat) != nsubs^2)
    #    stop("dimensions of outmat do not match nsubs and nseeds values")
    
    subsMap <- big.matrix(nvoxs-1, nsubs, type=type, shared=FALSE, ...)
    ALPHA <- 1/(nvoxs-2)
    voxs <- 1:nvoxs
    for (i in 1:nseeds) {
        .Call("CombineSubMapsMain", sub.cormaps, subsMap@address, as.double(i), 
              as.double(voxs[-cor_inds[i]]), as.double(nvoxs-1), as.double(nsubs))
        big_cor(x=subsMap, z=outmat, z_firstCol=dist_inds[i], z_lastCol=dist_inds[i])
        .Call("big_add_scalar_left", outmat, as.double(-1), as.double(dist_inds[i]), as.double(1));
    }
    rm(subsMap)
    gc(F)
    
    return(outmat)
}

compute_subdist_worker2_regress <- function(sub.cormaps, cor_inds, outmat, dist_inds,  
                                            design_mat, type="double", ...) 
{
    nsubs <- length(sub.cormaps)
    nvoxs <- ncol(sub.cormaps[[1]])
    nseeds <- nrow(sub.cormaps[[1]])
    if (nseeds != length(cor_inds) || nseeds != length(dist_inds))
        stop("length of inds doesn't match nrow of first sub.cormaps element")
    
    #if (is.null(outmat))
    #    outmat <- big.matrix(nsubs^2, nseeds, type=type, ...)
    #else if (ncol(outmat) != nseeds || nrow(outmat) != nsubs^2)
    #    stop("dimensions of outmat do not match nsubs and nseeds values")
    
    subsMap <- big.matrix(nsubs, nvoxs-1, type=type, shared=FALSE, ...)
    r_subsMap <- big.matrix(nsubs, nvoxs-1, type=type, shared=FALSE, ...)   # residuals
    voxs <- 1:nvoxs
    for (i in 1:nseeds) {
        .Call("CombineSubMapsTransSimpleMain", sub.cormaps, subsMap@address, as.double(i), 
              as.double(voxs[-cor_inds[i]]), as.double(nvoxs-1), as.double(nsubs))
        qlm_residuals(subsMap, design_mat, FALSE, r_subsMap)
        scale_fast(r_subsMap, to.copy=FALSE, byrows=TRUE)
        big_cor(x=r_subsMap, byrows=TRUE, z=outmat, 
                z_firstCol=dist_inds[i], z_lastCol=dist_inds[i])
        .Call("big_add_scalar_left", outmat, as.double(-1), as.double(dist_inds[i]), as.double(1));
    }
    
    rm(subsMap)
    gc(F)
    
    return(outmat)
}
