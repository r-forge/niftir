##-------------------
## Masking your niftiXd
##-------------------

#' @nord
read.mask <- function(fname, thresh=0) {
    x <- .Call("read_nifti", abspath(fname), 1, PACKAGE="niftir")
    ret <- as.vector(x$image)
    if (!is.null(thresh))
        ret>thresh
    else
        ret
}

#' @nord
write.mask <- function(...) write.nifti(...)


#' Access/set mask info
#' 
#' @name mask
#' @aliases mask<-
#' 
#' @usage
#'  mask(x)
#'  mask(x) <- y
#' 
#' @param x \code{niftiXd} or \code{big.niftiXd} object
#' 
#' @return logical vector for \code{mask(x)}
roxygen()

#' @nord
setGeneric('mask', function(x) standardGeneric('mask'))

#' @nord
setMethod('mask',
    signature(x='niftiXd'),
    function(x) x@mask
)

#' @nord
setMethod('mask',
    signature(x='big.niftiXd'),
    function(x) x@mask
)

#' @nord
setGeneric('mask<-', function(x, value) standardGeneric('mask<-'))

#' @nord
setMethod('mask<-',
    signature(x='niftiXd'),
    function(x, value) {
        if (!is.vector(value) && !is.logical(value))
            stop("mask must be vector and logical")
        x@mask <- value
    }
)

#' @nord
setMethod('mask<-',
    signature(x='big.niftiXd'),
    function(x, value) {
        if (!is.vector(value) && !is.logical(value))
            stop("mask must be vector and logical")
        x@mask <- value
    }
)

#' Prepare a mask to do masking
#' 
#' Intended as a function called by \code{do.mask}.
#' Ensures that input \code{new_mask} is a vector of logicals for masking.
#' 
#' If \code{orig_mask} is specified, then this implies that \code{new_mask}
#' will be masking an object later that has already been masked. The 
#' \code{new_mask} can then be either the same length as \code{orig_mask} or the 
#' length of the number of TRUE elements in \code{new_mask}.
#' 
#' @title Prepare mask
#' @author Zarrar Shehzad
#'
#' @usage prepare.mask(new_mask, orig_mask=NULL, thresh=0)
#'
#' @param new_mask vector or object that can be converted to a vector
#' @param orig_mask NULL or vector of logical elements
#' @param thresh if new_mask is numeric than threshold at this level
#'
#' @return vector of logicals
#'
#' @seealso \code{\link{do.mask}}
prepare.mask <- function(new_mask, orig_mask=NULL, thresh=0) {
    if (!is.numeric(new_mask) && !is.logical(new_mask))
        stop("new_mask must be a numeric or logical vector")
    if (!is.null(orig_mask) && !is.logical(orig_mask))
        stop("orig_mask must be a logical vector")
    
    ## enforce new_mask is a vector
    new_mask <- as.vector(new_mask)
    
    # 1. Ensure mask is logical
    if (!is.logical(new_mask))
        new_mask <- new_mask>thresh
    
    # 2. Check orig_mask
    if (is.null(orig_mask))
        return(new_mask)
    else if (!is.logical(orig_mask))
        stop("if orig mask specified, then must be vector of logical elements")
        
    # 3. Mask can have 2 possible lengths
    lenmask <- length(new_mask)
    ## a. = length of original mask
    if (lenmask == length(orig_mask)) {
        new_mask <- new_mask & orig_mask
    }    
    ## b. = length of total # of TRUE elements in orig_mask
    else if (lenmask == sum(orig_mask)) {
        tmp <- orig_mask
        tmp[tmp] <- new_mask
        new_mask <- tmp
    }
    ## error
    else {
        stop("mask must have length equal to length of mask slot in x or number 
            of columns of x")
    }
    
    return(new_mask)
}

#' Mask/unmask the voxels of your niftiXd object
#' 
#' Masking a niftiXd object involves including only those columns (voxels) that
#' are specified as TRUE or above a given threshold in the input \code{mask}.
#' 
#' The benefit of using this function over just \code{x[,mask]} is that it will
#' keep a record of the regions that were masked, allowing you to go back to the
#' original structure with \code{\link{do.unmask}} or easily save your file in 
#' the appropriate dimensions with \code{\link{write.nifti}}.
#' 
#' The \code{mask} input argument can have a length that is:
#' (1) equal to the total number of voxels in the original nifti image
#' (2) equal to the number of columns or elements in your input \code{x}
#' 
#' Note that columns or elements in a \code{niftiXd} correspond to different
#' voxels in 3D space.
#' 
#' The \code{do.unmask} function puts back the columns or vector elements that
#' were previously masked. The values of these new elements will be set to 0.
#' 
#' @name do.mask
#' @aliases do.unmask
#' @title Masking/unmasking your niftiXd object
#' @author Zarrar Shehzad
#' 
#' @usage
#'  do.mask(x, mask, thresh=0, output.prefix=NULL)
#'  do.unmask(x, return.niftiXd=TRUE)
#' 
#' @param x \code{niftiXd} object
#' @param mask A vector
#' @param thresh If the mask isn't logical, what values should we threshold it
#'  at (default: 0)
#' @param output.prefix if you have big.niftiXd object that is file-backed,
#'  you can specify this option
#' @param return.niftiXd For \code{do.unmask}, whether or not to return a 
#'  niftiXd type object (doesn't apply to big.niftiXd objects)
#' @param ... Additional options possible for \code{do.mask}
#' 
#' @return A masked \code{niftiXd} or \code{big.niftiXd} object
roxygen()

#' @nord
setGeneric('do.mask', function(x, mask, thresh=0, output.prefix=NULL, ...)
           standardGeneric('do.mask'))

#' @nord
setMethod('do.mask',
    signature(x='niftiXd', mask='character', output.prefix='missing'),
    function(x, mask, thresh=0) {
        mask <- read.nifti3d(mask)  ## assume 3D
        return(mask(x, mask, thresh))
    }
)

#' @nord
setMethod('do.mask',
    signature(x='big.niftiXd', mask='character'),
    function(x, mask, thresh=0, output.prefix=NULL, ...) {
        mask <- read.nifti3d(mask)  ## assume 3D
        return(mask(x, mask, thresh, output.prefix, ...))
    }
)

#' @nord
setMethod('do.mask',
    signature(x='big.nifti4d', mask='character'),
    function(x, mask, thresh=0, output.prefix=NULL, ...) {
        mask <- read.nifti3d(mask)  ## assume 3D
        return(mask(x, mask, thresh, output.prefix, ...))
    }
)

#' @nord
setMethod('do.mask',
    signature(x='big.nifti4d', mask='vector'),
    function(x, mask, thresh=0, output.prefix=NULL) {
        # Prepare mask
        mask <- prepare.mask(mask, x@mask, thresh)
        n <- sum(mask)
        
        # Check if file.backed
        if (is.filebacked(x) && is.null(output.prefix))
            warning("no output.prefix specified, will output a non-filebacked 
                matrix")
        
        # Get indices of mask
        inds <- which(mask[x@mask])
        
        # Set output stuff and create big matrix
        ylist <- niftir.big.matrix(nrow(x), n, output=output.prefix, type=typeof(x), separated=is.separated(x), shared=is.shared(x))
        
        # Call the function
        .Call("BigDeepCopyMain", x@address, ylist$bm@address, 
            as.double(1:nrow(x)), as.double(inds))
        
        # Create new class
        y <- new(class(x), mask=mask, backingfile=ylist$bf, descriptorfile=ylist$df)
        y@address <- ylist$bm@address
        
        # Add any other slot names that might have been missed
        slot_names <- setdiff(slotNames(x), 
                                c("mask", "address", "backingfile", "descriptorfile"))
        for (slot_name in slot_names)
            slot(y, slot_name) <- slot(x, slot_name)
        
        return(y)
    }
)

#' @nord
setMethod('do.mask',
    signature(x='niftiXd', mask="vector", output.prefix="missing"),
    function(x, mask, thresh=0) {  
        mask <- prepare.mask(mask, x@mask, thresh)
        x@mask <- mask
        x@.Data <- x[,mask]
        return(x)
    }
)

#' @nord
setGeneric('do.unmask', function(x, return.niftiXd=TRUE)
           standardGeneric('do.unmask'))

#' @nord
setMethod('do.unmask',
    signature(x='big.niftiXd', return.niftiXd='missing'),
    function(x) {
        stop("unmasking is not yet supported with a big.niftiXd object")
})

#' @nord
setMethod('do.unmask',
    signature(x='nifti3d'),
    function(x, return.niftiXd=TRUE) {
        # Check if need to unmask
        lenmask <- length(x@mask)
        if (lenmask == length(x)) {
            if (return.niftiXd)
                return(x)
            else
                return(x@.Data)
        }
        
        x@header$dim <- x@header$dim[1:3]
        newVec <- vector("numeric", prod(x@header$dim))
        newVec[x@mask] <- x@.Data
        
        if (return.niftiXd) {
            x@.Data <- newVec
            x@mask <- rep(TRUE, lenmask)
            return(x)
        } else {
            return(newVec)
        }
})

#' @nord
setMethod('do.unmask',
    signature(x='nifti4d'),
    function(x, return.niftiXd=TRUE) {
        # Check if need to unmask
        lenmask <- length(x@mask)
        if (lenmask == ncol(x)) {
            if (return.niftiXd)
                return(x)
            else
                return(x@.Data)
        }
        
        x@header$dim <- x@header$dim[1:4]
        newMat <- matrix(0, x@header$dim[4], prod(x@header$dim[1:3]))
        newMat[,x@mask] <- x@.Data
        
        if (return.niftiXd) {
            x@.Data <- newMat
            x@mask <- rep(TRUE, lenmask)
            return(x)
        } else {
            return(newMat)
        }
})


#-------------------
# Indices of your niftiXd
#-------------------

#' Get the indices of your mask
#' 
#' @name indices
#' @aliases indices-methods
#' @author Zarrar Shehzad
#' 
#' @param x \code{niftiXd} object or mask
#' @param i mask that might be applied to \code{x]}
#'
#' @return vector of column or vector indices in \code{x}
roxygen()

#' @nord
setGeneric('indices', function(x, i) standardGeneric('indices'))

#' @nord
setMethod('indices',
    signature(x='logical', i="logical"),
    function(x, i) {
        i <- prepare.mask(i, x)
        i <- which(i[x])
        indices(x, i)
})

#' @nord
setMethod('indices',
    signature(x='logical', i="integer"),
    function(x, i) {
        if (min(i)<1 || max(i)>length(x))
            stop("Illegal index usage")
        return(i)
})

#' @nord
setMethod('indices',
    signature(x='niftiXd', i="logical"),
    function(x, i) {
        indices(x@mask, i)
})

#' @nord
setMethod('indices',
    signature(x='big.niftiXd', i="logical"),
    function(x, i) {
        indices(x@mask, i)
})
