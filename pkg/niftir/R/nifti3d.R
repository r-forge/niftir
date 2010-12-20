#----------------------
# as.nifti3d FUNCTIONS
#----------------------

#' Create a nifti3d object
#'
#' @name as.nifti3d-methods
#'
#' @aliases as.nifti3d,nifti-method
#'  as.nifti3d,array,list-method
#'  as.nifti3d,vector,list,logical-method
#'
#' @sectionMethods
#'  \describe{
#'      \item{\code{signature(x="nifti")}}{...}
#'      \item{\code{signature(x="array", header="list")}}{...}
#'      \item{\code{signature(x="vector", header="list", mask="logical")}}{...}
#'  }
#'
#' @seealso \code{\link{as.nifti3d}}
#'
#' @keywords methods
roxygen()

#' Create a nifti3d object
#'
#' @name as.nifti3d
#'
#' @usage
#'  as.nifti3d(x)   # when x is a nifti object
#'  as.nifti3d(x, header)   # when x is an array
#'  as.nifti3d(x, header, mask) #when x is a vector
#'
#' @author Zarrar Shehzad
#' 
#' @param x 3D \code{nifti}, 3D \code{array}, or \code{vector}
#' @param header list of header attributes (required when \code{x} is an
#'  \code{array} or \code{vector})
#' @param mask logical vector specifying which voxels from 3D image are
#'  specified with \code{x} (only required when \code{x} is a \code{vector})
#'
#' @return \code{nifti4d} object
#' 
#' @seealso \code{\link{as.nifti}}, \code{\link{as.nifti3d-methods}}
#'
#' @examples
#'  as.nifti3d(array(0, c(10,10,10)), create.header())
#'  as.nifti3d(nifti(0, dim=c(10,10,10)))   # should give same thing as above
#' 
#' @keywords methods
roxygen()

#' @nord
setGeneric('as.nifti3d', 
    function(x, header, mask) standardGeneric('as.nifti3d')
)

#' @nord
setMethod('as.nifti3d',
    signature(x='nifti', header='missing', mask='missing'),
    function(x) as.nifti3d(x@.Data, header(x))
)

#' @nord
setMethod('as.nifti3d',
    signature(x='array', header='list', mask='missing'),
    function(x, header) {
        if (length(dim(x)) != 3)
            stop("dimensions not equal to 3")
        
        # 1. Create 1D vector
        header$dim <- dim(x)  # ensure header/image consistency
        x <- as.vector(x)
        
        # 2. Create new class
        vec <- new("nifti3d", header=header, mask=rep(TRUE, length(x)))
        vec@.Data <- x
        
        # 3. Set names with 'x:y:z' for each element
        coords <- coords(vec)
        rn <- apply(coords, 1, paste, collapse=":")
        names(vec) <- rn
        
        gc()
        
        return(vec)
    }
)

#' @nord
setMethod('as.nifti3d',
    signature(x='vector', header='list', mask='logical'),
    function(x, header, mask) {
        # 1. Check input
        if (length(header$dim) != 3)
            stop("dimensions of header attribute not equal to 3")
        if (prod(header$dim) != length(mask))
            stop("mask must have length equal to total number of voxels in
                 3D nifti image")
        if (length(x) != sum(mask))
            stop("mask must have as many TRUE elements as length of input x")
        
        # 2. Create new class
        vec <- new("nifti3d", header=header, mask=mask)
        vec@.Data <- x
        
        # 3. Set names with 'x:y:z' for each element
        coords <- coords(vec)
        rn <- apply(coords, 1, paste, collapse=":")
        names(vec) <- rn
        
        gc()
        
        return(vec)
    }
)



#' @nord
setGeneric('is.nifti3d', 
    function(x) standardGeneric('is.nifti3d')
)

#' @nord
setMethod('is.nifti3d',
    signature(x='nifti3d'),
    function(x) return(TRUE)
)

#' @nord
setMethod('is.nifti3d',
    definition=function(x) return(FALSE)
)


#----------------------
# read.nifti3d FUNCTION
#----------------------

#' Read in a nifti3d object from a file
#'
#' @usage read.nifti3d(fname)
#'
#' @author Zarrar Shehzad
#' 
#' @param fname character specifying path to analyze/nifti file
#'
#' @return \code{nifti4d} object
#' 
#' @seealso \code{\link{read.nifti}}, \code{\link{as.nifti3d}}
#'
#' @examples
#'  fname <- file.path(system.file("data", package="niftir"), "test.nii.gz") # 2x2x2 size file
#'  # or fname <- file.choose()
#'  nim <- read.nifti3d(fname)
#' 
#' @keywords methods
read.nifti3d <- function(fname) {
    x <- read.nifti(fname)
    as.nifti3d(x)
}

