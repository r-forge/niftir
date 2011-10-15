# Higher-level wrappers around certain R functions

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
    vec <- reho(bmat, verbose=verbose, parallel=parallel, ...)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        write.nifti(vec, hdr, mask, outfile=out_file, overwrite=overwrite)
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(vec)
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
    bmats <- load_and_mask_func_data(func_files, mask, type="double", shared=shared)
    
    vcat(verbose, "...kendalling")
    verbosity <- ifelse(verbose, 2, 0)
    vec <- kendall(bmats, blocksize, verbosity=verbosity, parallel=parallel, ...)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        write.nifti(vec, hdr, mask, outfile=out_file, overwrite=overwrite, odt="float")
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(vec)
    }
}

