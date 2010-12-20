#----------------------
# as.nifti3d FUNCTIONS
#----------------------

#' Create a nifti4d object
#'
#' @name as.nifti4d-methods
#'
#' @aliases as.nifti4d,nifti-method
#'  as.nifti4d,big.nifti4d-method
#'  as.nifti4d,array,list-method
#'  as.nifti4d,matrix,list,logical-method
#'
#' @sectionMethods
#'  \describe{
#'      \item{\code{signature(x="nifti")}}{...}
#'      \item{\code{signature(x="big.nifti4d")}}{...}
#'      \item{\code{signature(x="array", header="list")}}{...}
#'      \item{\code{signature(x="matrix", header="list", mask="logical")}}{...}
#'  }
#'
#' @seealso \code{\link{as.nifti4d}}
#'
#' @keywords methods
roxygen()

#' Create a nifti4d object
#'
#' @name as.nifti4d
#'
#' @usage
#'  as.nifti4d(x)   # when x is a nifti object or big.nifti4d object
#'  as.nifti4d(x, header)   # when x is an array
#'  as.nifti4d(x, header, mask) # when x is a matrix
#'
#' @author Zarrar Shehzad
#' 
#' @param x 4d \code{nifti}, 4d \code{array}, \code{big.nifti4d}, or \code{matrix}
#' @param header list of header attributes (required when \code{x} is an
#'  \code{array} or \code{matrix})
#' @param mask logical vector specifying which voxels from 3D image are
#'  specified with \code{x} (only required when \code{x} is a \code{matrix})
#'
#' @return \code{nifti4d} object
#' 
#' @seealso \code{\link{as.nifti}}, \code{\link{as.nifti4d-methods}}
#'
#' @examples
#'  as.nifti4d(array(0, c(10,10,10,10)), create.header())
#'  as.nifti4d(nifti(0, dim=c(10,10,10,10)))   # should give same thing as above
#' 
#' @keywords methods
roxygen()

#' @nord
setGeneric('as.nifti4d', 
    function(x, header, mask) standardGeneric('as.nifti4d')
)

#' @nord
setMethod('as.nifti4d',
    signature(x='nifti', header='missing', mask='missing'),
    function(x) as.nifti4d(x@.Data, x@header)
)

#' @nord
setMethod('as.nifti4d',
    signature(x='array', header='list', mask='missing'),
    function(x, header) {
        if (length(dim(x)) != 4)
            stop("dimensions not equal to 4")
        
        # 1. Create 2D matrix
        header$dim <- dim(x)  # ensure header/image consistency
        dim(x) <- c(prod(dim(x)[1:3]), dim(x)[4])  # rows = voxels & cols = timepoints
        x <- aperm(x)   # rows = timepoints & cols = voxels
        
        # 2. Create new class
        mat <- new("nifti4d", header=create.header(header), mask=rep(TRUE, ncol(x)))
        mat@.Data <- x

        # 3. Set dimnames with 'x:y:z' for each column name
        coords <- coords(mat)
        names(dim(mat)) <- c("voxels", "timepoints")
        rn <- apply(coords, 1, paste, collapse=":")
        dimnames(mat) <- list(timepoints=1:nrow(x), voxels=rn)
        
        gc()
        
        return(mat)
    }
)

#' @nord
setMethod('as.nifti4d',
    signature(x='matrix', header='list', mask='logical'),
    function(x, header, mask) {
        # 1. Check input
        if (length(header$dim) != 4)
            stop("dimensions of header attribute not equal to 4")
        if (prod(header$dim[1:3]) != length(mask))
            stop("mask must have length equal to total number of voxels in
                 nifti image")
        if (ncol(x) != sum(mask))
            stop("mask must have as many TRUE elements as columns in input x")
        
        # 2. Create new class
        mat <- new("nifti4d", header=create.header(header), mask=mask)
        mat@.Data <- x
        
        # 3. Set dimnames with 'x:y:z' for each column name
        coords <- coords(mat)
        names(dim(mat)) <- c("voxels", "timepoints")
        rn <- apply(coords, 1, paste, collapse=":")
        dimnames(mat) <- list(timepoints=1:nrow(x), voxels=rn)
        
        gc()
        
        return(mat)
    }
)

#' @nord
setMethod('as.nifti4d',
    signature(x='big.nifti4d'),
    function(x) as.nifti(x[,], x@header, x@mask)
)


#' @nord
setGeneric('is.nifti4d', 
    function(x) standardGeneric('is.nifti4d')
)

#' @nord
setMethod('is.nifti4d',
    signature(x='nifti4d'),
    function(x) return(TRUE)
)

#' @nord
setMethod('is.nifti4d',
    definition=function(x) return(FALSE)
)



#----------------------
# read.nifti4d FUNCTION
#----------------------

#' Read in a nifti4d object from a file
#'
#' @usage read.nifti4d(fname)
#'
#' @author Zarrar Shehzad
#' 
#' @param fname character specifying path to analyze/nifti file
#'
#' @return \code{nifti4d} object
#' 
#' @seealso \code{\link{read.nifti}}, \code{\link{as.nifti4d}}
#'
#' @examples
#'  # TODO
#' 
#' @keywords methods
read.nifti4d <- function(fname) {
    x <- read.nifti(fname)
    as.nifti4d(x)
}

#' @nord
setGeneric('voxapply', 
    function(x, FUN, header, mask, ...) standardGeneric('voxapply')
)

#' @nord
setMethod('voxapply',
    signature(x='matrix', FUN="function", header='list', mask='logical'),
    function(x, FUN, header, mask, ...) {
        vec <- apply(x, 2, FUN, ...)
        header$dim <- header$dim[1:3]
        header$pixdim <- header$pixdim[1:3]
        as.nifti3d(vec, header, mask)
    }
)

#' @nord
setMethod('voxapply',
    signature(x='nifti4d', FUN="function", header='missing', mask='missing'),
    function(x, FUN, ...) voxapply(x@.Data, FUN, x@header, x@mask)
)
