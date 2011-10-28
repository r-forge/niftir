# Higher-level wrappers around certain R functions

roi_mean_wrapper <- function(func_file, roi_file, mask_file=NULL, 
                              out_file=NULL, outtype="nifti", 
                              to_return=FALSE, overwrite=FALSE, 
                              verbose=TRUE)
{
    vcat(verbose, "Averaging mean signal in '%s' with '%s' ROIs", 
            basename(func_file), basename(roi_file))
    
    if (!is.character(func_file) || !file.exists(func_file))
        stop("Could not find functional file ", func_file)
    if (!is.character(roi_file) || !file.exists(roi_file))
        stop("Could not find roi file ", roi_file)
    if (!is.null(out_file) && file.exists(out_file) && !overwrite) {
        vcat(verbose, "File '%s' already exists, not re-running", out_file)
        return(NULL)
    }
    if (is.null(out_file) && !to_return)
        stop("You haven't specified an output file or asked to return the data")
    if (!(outtype %in% c("nifti", "text")))
        stop("output type can only be nifti or text")
    
    vcat(verbose, "...reading data")
    func <- read.big.nifti4d(func_file)
    rois <- read.mask(roi_file, NULL)
    hdr <- read.nifti.header(func_file)
    if (!is.null(mask_file)) {
        mask <- read.mask(mask_file) & rois!=0
    } else {
        mask <- rois!=0
    }
    
    vcat(verbose, "...masking")
    func <- do.mask(func, mask)
    rois <- rois[mask]
    
    vcat(verbose, "...averaging")
    new_func <- roi_mean(func, rois)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        if (outtype == "nifti") {
            hdr$dim <- dim(new_func)
            hdr$pixdim <- c(new_func[4], 1)
            write.nifti(new_func, hdr, odt="float", 
                        outfile=out_file, overwrite=overwrite)
        } else if (outtype == "text") {
            write.table(new_func, file=out_file, quote=FALSE, 
                        row.names=F, col.names=unique(rois))
        }
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(new_func)
    }
}

wrap_gcor <- function(func_file, mask_file, out_file=NULL, 
                      blocksize=0, memlimit=4, 
                      to_return=FALSE, overwrite=FALSE, 
                      verbose=TRUE, parallel=FALSE, shared=parallel, 
                        ...) 
{
    vcat(verbose, "Running global correlation for '%s'", func_file)
    
    if (!is.character(func_file) || !file.exists(func_file))
        stop("Could not find functional file ", func_file)
    if (!is.character(mask_file) || !file.exists(mask_file))
        stop("Could not find mask file ", mask_file)
    if (!is.null(out_file) && file.exists(out_file) && !overwrite) {
        vcat(verbose, "File '%s' already exists, not re-running", out_file)
        return(NULL)
    }
    
    vcat(verbose, "...reading mask")
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    
    vcat(verbose, "...calculating memory demands")
    tmp <- read.nifti.header(func_file)
    if (length(tmp$dim) != 4)
        vstop("functional file '%s' must be 4D", func_file)
    blocksize <- get_gcor_limit(blocksize, memlimit, sum(mask), tmp$dim[[4]], verbose)
    
    vcat(verbose, "...reading functional")
    bmat <- load_and_mask_func_data(func_file, mask, type="double", shared=shared)
    
    vcat(verbose, "...correlating")
    vec <- gcor(bmat, blocksize, verbosity=verbose*1, parallel=parallel, shared=shared, 
                ...)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        write.nifti(vec, hdr, mask, outfile=out_file, overwrite=overwrite)
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(vec)
    }
}


wrap_reho <- function(func_file, mask_file, out_file=NULL, 
                      to_return=FALSE, overwrite=FALSE, 
                      verbose=TRUE, parallel=FALSE, shared=parallel, 
                      ...) 
{
    vcat(verbose, "Running regional homogeneity for '%s'", func_file)
    
    if (!is.character(func_file) || !file.exists(func_file))
        stop("Could not find functional file ", func_file)
    if (!is.character(mask_file) || !file.exists(mask_file))
        stop("Could not find mask file ", mask_file)
    if (!is.null(out_file) && file.exists(out_file) && !overwrite) {
        vcat(verbose, "File '%s' already exists, not re-running", out_file)
        return(NULL)
    }
    
    vcat(verbose, "...reading mask")
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    
    vcat(verbose, "...reading functional")
    bmat <- load_and_mask_func_data3(func_file, mask, type="double", shared=shared)
    
    vcat(verbose, "...rehoing")
    res <- reho(bmat, verbose=verbose, parallel=parallel, ...)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        write.nifti(res$reho, res$hdr, res$mask, outfile=out_file, overwrite=overwrite)
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(res)
    }
}


wrap_kendall <- function(func_files, mask_file, out_file=NULL, to_return=FALSE, 
                            overwrite=FALSE, verbose=TRUE, parallel=FALSE, 
                            shared=parallel, memlimit=4, ...)
{
    vcat(verbose, "Running kendall's W")
    progress <- ifelse(verbose, "text", "none")
    
    if (!is.character(func_files) && length(func_files) < 2)
        stop("func_files must be a vector of at least 2 filenames")
    for (func_file in func_files) {
        if (!is.character(func_file) || !file.exists(func_file))
            stop("Could not find functional file ", func_file)
    }
    if (!is.character(mask_file) || !file.exists(mask_file))
        stop("Could not find mask file ", mask_file)    
    if (!is.null(out_file) && file.exists(out_file) && !overwrite) {
        vcat(verbose, "File '%s' already exists, not re-running", out_file)
        return(NULL)
    }
    
    vcat(verbose, "...reading mask")
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    
    vcat(verbose, "...getting # of time-points for functional data")
    ntpts <- laply(func_files, function(x) {
        hdr <- read.nifti.header(x)
        if (length(hdr$dim) != 4) {
            vstop("Input file '%s' must be 4 dimensions but is %i dimensional", 
                  x, length(hdr$dim))
        }
        return(hdr$dim[[4]])
    }, .progress=progress)
    
    vcat(verbose, "...calculating memory demands")
    blocksize <- get_kendall_limit(0, memlimit, sum(mask), ntpts, verbose)
    
    vcat(verbose, "...reading %i functionals", length(func_files))
    reader <- gen_big_reader("nifti4d", type="double", shared=shared)
    funclist <- load_and_mask_func_data2(func_files, reader, mask=mask, 
                                         verbose=verbose,  
                                         type="double", shared=shared)
    
    vcat(verbose, "...checking functionals")
    checks <- check_func_data(func_files, funclist, verbose=verbose, parallel=parallel)
    
    vcat(verbose, "...kendalling")
    verbosity <- ifelse(verbose, 2, 0)
    vec <- kendall(funclist, blocksize, verbosity=verbosity, parallel=parallel, ...)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        write.nifti(vec, hdr, mask, outfile=out_file, overwrite=overwrite, odt="float")
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(vec)
    }
}

wrap_glm <- function(func_files, mask_file, ev_file, contrast_file, 
                     outdir, overwrite=FALSE, 
                     blocksize=0, memlimit=4, 
                     verbose=TRUE, parallel=FALSE, shared=parallel, 
                     ztransform=FALSE) 
{
    vcat(verbose, "Running GLM")
    
    vcat(verbose, "...setup and checks")
    progress <- ifelse(verbose, "text", "none")
    if (!is.character(func_files) && length(func_files) < 2)
        stop("func_files must be a vector of at least 2 filenames")
    for (func_file in func_files) {
        if (!is.character(func_file) || !file.exists(func_file))
            stop("Could not find functional file ", func_file)
    }
    if (!is.character(mask_file) || !file.exists(mask_file))
        stop("Could not find mask file ", mask_file)   
    if (!is.character(ev_file) || !file.exists(ev_file)) 
        stop("Could not find EV file ", ev_file)
    if (!is.character(contrast_file) || !file.exists(contrast_file)) 
        stop("Could not find contrast file ", contrast_file)
    
    vcat(verbose, "...reading and setting up EV and contrast files")
    evs <- as.matrix(read.table(ev_file, header=T))
    cons <- as.matrix(read.table(contrast_file, header=T))
    contrast_names <- rownames(cons)
    if (is.numeric(contrast_names))
        stop("must have row names for contrast matrix")
    if (ncol(evs) != ncol(cons)) {
        vstop(paste("# of columns in EV file '%s' not the same as # of columns",  
                    "in contrast file '%s'"), ev_file, contrast_file)
    }
    tmp <- big.matrix(nrow(evs), ncol(evs), type="double", shared=shared)
    tmp[,] <- evs; evs <- tmp; rm(tmp)
    #tmp <- big.matrix(nrow(cons), ncol(cons), type="double", shared=shared)
    #tmp[,] <- cons; cons <- tmp; rm(tmp)
    k <- qlm_rank(evs)
    if (k < ncol(evs))
        vstop("EV file '%s' is rank deficient", ev_file)
    nevs <- ncol(evs)
    ncons <- nrow(cons)
    
    vcat(verbose, "...reading mask")
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    nvoxs <- sum(mask)
    
    vcat(verbose, "...setting up output")
    outdir <- abspath(outdir)
    if (!is.null(outdir) && file.exists(outdir)) {
        if (!overwrite) {
            vcat(verbose, "Directory '%s' already exists, not re-running", outdir)
            return(NULL)
        } else {
            vcat(verbose, "Removing directory '%s'", outdir)
            system(sprintf("rm %s/infuncs/*", outdir))
            system(sprintf("rmdir %s/infuncs", outdir))
            system(sprintf("rm %s/*", outdir))
            system(sprintf("rmdir %s", outdir))
        }
    } else {
        dir.create(outdir)
    }
    ## copy evs
    file.copy(ev_file, file.path(outdir, "model_evs.txt"))
    ## copy contrasts
    file.copy(contrast_file, file.path(outdir, "model_contrasts.txt"))
    ## copy mask
    file.copy(mask_file, file.path(outdir, "mask.nii.gz"))
    ## soft-link functionals
    dir.create(file.path(outdir, "infuncs"))
    for (i in 1:length(func_files)) {
        file.symlink(func_files[i], file.path(outdir, "infuncs", 
                                           sprintf("func%04i.nii.gz", i)))
    }
    ## create output matrices
    tmats <- lapply(1:ncons, function(i) {
        bfile <- sprintf("tvals_%02i.bin", i)
        dfile <- sprintf("tvals_%02i.desc", i)
        big.matrix(nvoxs, nvoxs, type="double", 
                   backingfile=bfile, descriptorfile=dfile, 
                   backingpath=outdir)
    })
    
    vcat(verbose, "...getting # of time-points for functional data")
    ntpts <- laply(func_files, function(x) {
        hdr <- read.nifti.header(x)
        if (length(hdr$dim) != 4) {
            vstop("Input file '%s' must be 4 dimensions but is %i dimensional", 
                  x, length(hdr$dim))
        }
        return(hdr$dim[[4]])
    }, .progress=progress)
    
    vcat(verbose, "...calculating memory demands")
    blocksize <- get_glm_limit(0, memlimit, nvoxs, ntpts, 
                               nevs, ncons, verbose)
    
    vcat(verbose, "...reading %i functionals", length(func_files))
    reader <- gen_big_reader("nifti4d", type="double", shared=shared)
    funclist <- load_and_mask_func_data2(func_files, reader, mask=mask, 
                                         verbose=verbose,  
                                         type="double", shared=shared)
    
    vcat(verbose, "...checking functionals")
    checks <- check_func_data(func_files, funclist, verbose=verbose, parallel=parallel)
    
    # loop through a set of voxels
    ## calculate connectivity maps
    ## create temporary contrasts
    ## get fit for given seed map
    ## get contrasts (results)
    ## save contrasts and clear file-backed stuff
    vcat(verbose, "...GLMing")
    start.time <- Sys.time()
    vox_glm(funclist, evs, cons, blocksize, outmats=tmats, bp=outdir, 
            verbose=verbose, parallel=parallel, shared=shared, 
            ztransform=ztransform)
    end.time <- Sys.time()
    vcat(verbose, "GLMing is done! It took: %.2f minutes\n", 
         as.numeric(end.time-start.time, units="mins"))
    
    invisible(tmats)
}
