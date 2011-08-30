# Higher-level wrappers around certain R functions

wrap_gcor <- function(func_file, mask_file, out_file=NULL, to_return=FALSE, 
                        overwrite=FALSE, verbose=TRUE, ...) 
{
    vcat(verbose, "Running global correlation for '%s'", func_file)
    
    if (!is.character(func_file) || !file.exists(func_file))
        stop("Could not find functional file ", func_file)
    if (!is.character(mask_file) || !file.exists(mask_file))
        stop("Could not find mask file ", mask_file)
    
    vcat(verbose, "...reading mask")
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    
    vcat(verbose, "...reading functional")
    bmat <- load_and_mask_func_data(func_file, mask)
    
    vcat(verbose, "...correlating")
    vec <- gcor(bmat, verbose=verbose, ...)
    
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
                            overwrite=FALSE, verbose=TRUE, ...)
{
    vcat(verbose, "Running kendall's W")
    
    if (!is.list(func_files))
        stop("func_files must be a list")
    for (func_file in func_files) {
        if (!is.character(func_file) || !file.exists(func_file))
            stop("Could not find functional file ", func_file)
    }
    if (!is.character(mask_file) || !file.exists(mask_file))
        stop("Could not find mask file ", mask_file)
    
    vcat(verbose, "...reading mask")
    mask <- read.mask(mask_file)
    hdr <- read.nifti.header(mask_file)
    
    vcat(verbose, "...reading %i functionals", length(func_files))
    bmats <- load_and_mask_func_data(func_files, mask)
    
    vcat(verbose, "...kendalling")
    vec <- kendall(bmats, ...)
    
    if (!is.null(out_file)) {
        vcat(verbose, "...writing file '%s'", out_file)
        write.nifti(vec, hdr, mask, outfile=out_file, overwrite=overwrite)
    }
    
    if (to_return) {
        vcat(verbose, "...returning data")
        return(vec)
    }
}

