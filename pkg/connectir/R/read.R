# Wrapper read functions
read.big.tab <- function(file, ...) read.big.matrix(file, sep='\t', ...)

read.big.space <- function(file, ...) read.big.matrix(file, sep=' ', ...)

read.big.csv <- function(file, ...) read.big.matrix(file, sep=',', ...)

read.big.matlab <- function(file, ...) {
    library(R.matlab)
    m <- readMat(file)
    if (length(m) > 1)
        stop(sprintf("matlab input '%s' can't have more than 1 variable, but %i found", 
                        fname, length(m)))
    as.big.matrix(m[[1]], ...)
}

# Read input data as big matrix
gen_big_reader <- function(intype, ...) {
    choices <- c("nifti4d", "tab", "space", "csv", "matlab")
    if (!(intype %in% choices))
        stop("unrecognized input type: ", intype)
    
    args <- list(...)
    fun <- function(x, ...) {
        # Checks
        if (is.big.matrix(x))
            return(x)
        else if (!is.character(x))
            stop("input must be a character or big.matrix and not ", class(x))
        
        # Read
        args$file <- x
        mat <- do.call(sprintf("read.big.%s", intype), args)
        
        gc(FALSE, TRUE)
        return(mat)
    }
    
    return(fun)
}

# Detect the type of files
detect_ftypes <- function(fnames, force.type=NULL, verbose=TRUE) {
    df <- data.frame(
        extensions = c(".nii.gz", ".nii", ".hdr", ".img", ".mat", ".txt", ".tab", ".csv"), 
        formats = c("nifti4d", "nifti4d", "nifti4d", "nifti4d", "matlab", "space", "tab", "csv")
    )
    
    if (!is.null(force.type) && !(force.type %in% df$formats))
        stop("unrecognized file type ", force.type)
    
    formats <- sapply(fnames, function(fname) {
        find <- sapply(df$extensions, function(ext) grepl(ext, fname))
        if (any(find))
            return(as.character(df$formats[find]))
        else if (is.null(force.type))
            stop("unrecognized extension in ", fname)
    })
    
    if (any(formats != formats[1]))
        warning("not all extensions are ", formats[1], " like in ", fnames[1])
    if (!is.null(force.type) && any(formats != force.type))
        vcat(verbose, "not all extensions are the same as the forced one ", force.type)
    
    ftype <- ifelse(is.null(force.type), force.type, formats[1])
    
    return(ftype)
}

# Automatically determine the appropriate mask for data
## computes variance @ each voxel
## if var=0, then vox=0, otherwise vox=1
automask <- function(x, cols=NULL, na.rm=FALSE) {
    vs <- colvar(x, cols, na.rm)
    return(vs!=0)
}

# Create common mask across participants
## '...' => for as.big.matrix creation
## exclude.thresh => 
overlap_automasks <- function(xs, read_fun, verbose=FALSE, parallel=FALSE, na.rm=FALSE, 
                              exclude.thresh=0, ...) 
{
    if (!is.list(xs) && !is.vector(xs))
        stop("input 'xs' must be a vector or list")
    if (!is.function(read_fun))
        stop("input 'read_fun' must be a function")
    progress <- ifelse(verbose, "text", "none")
    n <- length(xs)
    
    masks <- laply(xs, function(x) {
        x <- read_fun(x, shared=parallel, ...)
        m <- automask(x, na.rm=na.rm)*1
        rm(x); gc(FALSE, TRUE)
        return(m)
    }, .progress=progress, .parallel=parallel)
    
    nas <- is.na(masks)
    if (any(nas)) {
        vcat(verbose, "%i NaNs found...setting to 0", sum(nas))
        masks[nas] <- 0
    }
    
    nnodes <- ncol(masks)
    sub.ns <- vector("numeric", n)
    subs.nbad <- alply(masks, 2, function(x) which(x!=1))
    for (i in 1:nnodes) 
        sub.ns[subs.nbad[[i]]] <- sub.ns[subs.nbad[[i]]] + 1
    exclude.subs <- which(sub.ns/nnodes > exclude.thresh)
    
    overlap <- colMeans(masks[-exclude.subs,])
    
    mask <- overlap == 1
    vcat(verbose, "mask has %i good nodes and %i bad ones", sum(mask), sum(!mask))
    
    return(list(mask=mask, overlap=overlap, sub.masks=masks))
}

# Load data
load_and_mask_func_data2 <- function(xs, read_fun, mask=NULL, verbose=FALSE, ...)
{
    if (!is.list(xs) && !is.vector(xs))
        stop("input 'xs' must be a vector or list")
    if (!is.function(read_fun))
        stop("input 'read_fun' must be a function")
    progress <- ifelse(verbose, "text", "none")
    
    if (!is.null(mask) && is.logical(mask))
        mask <- which(mask)
    
    vcat(verbose, "reading data")
    dat.list <- llply(xs, function(x) {
        x <- read_fun(x, ...)
        if (!is.null(mask)) {
            z <- deepcopy(x, cols=mask, ...)
            rm(x); gc(FALSE, TRUE)
            x <- z
            rm(z); gc(FALSE, TRUE)
        }
        y <- scale(x, to.copy=TRUE, ...)
        rm(x); gc(FALSE, TRUE)
        return(y)
    }, .progress=progress, .parallel=FALSE)
    
    gc(FALSE, TRUE)
    
    return(dat.list)
}

check_func_data <- function(dat.list, verbose=FALSE, parallel=FALSE) 
{
    vcat(verbose, "checking data")
    nc <- ncol(dat.list[[1]])
    n <- length(dat.list)
    rets <- llply(1:n, function(i) {
        dat <- dat.list[[i]]
        fname <- xs[[i]]
        
        # Everything must have same # of voxels
        if (nc != ncol(dat)) {
            vcat(verbose, 
                 "%s must have the same # of nodes (%i vs %i) as other datasets", 
                  fname, nc, ncol(dat))
            return(1)
        }
        
        # There can't be any NaNs
        col.nas <- colna(dat.list[[i]])>0
        if (any(col.nas)) {
            w <- paste(which(col.nas), collapse=", ")
            vcat(verbose, "%s has NaNs in nodes %s", fname, w)
            return(2)
        }
        
        return(0)
    }, .progress=progress, .parallel=parallel)
    
    unlist(rets)
}
