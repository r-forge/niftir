# TODO: figure out more decriptive class names for the big.nifti5d and big.nifti6d
# TODO: check if performance decrements from level of recursion in class defs

#' Class containing a list of header attributes
#' 
#' @title header class
#' @slot header A list that holds attributes of a nifti image
setClass(
    Class = "header", 
    representation = representation(
        header = "list"
    )
)


#' Class returned when reading a nifti/analyze file or from \code{\link{nifti}}/\code{\link{as.nifti}}
#' 
#' @title nifti/analyze class
#' @slot header A list that holds attributes of a nifti image
#' @slot .Data This holds the image data
setClass(
    Class = "nifti", 
    representation = representation(),
    contains = c("header", "array")
)


#' Not really used on its own (but the class has friends that extend it)
#' 
#' @slot mask A logical vector (represents the voxels in the original 3D image whose values are given in \code{.Data} or \code{address}); its length must equal the product of the dim attribute in the header
setClass(
    Class = "niftiXd",
    representation = representation(
        mask = "vector"
    )
)


#' Represents a 3D nifti image as a 1D vector
#' 
#' @title 3d nifti => 1d vector
#' @slot header A list that holds attributes of a nifti image
#' @slot mask A logical vector (represents the voxels in the original 3D image whose values are given in \code{.Data}); its length must equal the product of the dim attribute in the header
#' @slot .Data This holds the image data in 1D
setClass(
    Class = "nifti3d",
    representation = representation(),
    contains = c("header", "niftiXd", "vector")
)


#' Represents a 4D nifti image as a 2D matrix
#' 
#' Assuming that the 4D image is akin to 3D+time, then the columns represent the voxels in 3D while the rows represent the different timepoints.
#' 
#' @title 4d nifti => 2d matrix
#' @slot header A list that holds attributes of a nifti image
#' @slot mask A logical vector (represents the voxels in the original 4D image whose values are given in each column of \code{.Data}); its length must equal the product of the dim attribute in the header
#' @slot .Data This holds the image data in 2D
setClass(
    Class = "nifti4d",
    representation = representation(),
    contains = c("header", "niftiXd", "matrix")
)


#' Not really used on its own (but the class has friends that extend it)
#' 
#' @slot header A list that holds attributes of a nifti image
#' @slot mask A logical vector (represents the voxels in the original 3D image whose values are given in \code{.Data} or \code{address}); its length must equal the product of the dim attribute in the header
#' @slot address Object of class \code{externalptr} points to the memory location of the C++ data structure.
#' @slot bpath backing path if filebacked big matrix
#' @slot dfile descriptor filename
setClass(
    Class = "big.niftiXd",
    representation = representation(
        mask = "vector",
        backingfile = "character",
        descriptorfile = "character"
    ),
    contains = c("big.matrix")
)


#' Represents a 4D nifti image as a 2D big matrix
#' 
#' Assuming that the 4D image is akin to 3D+time, then the columns represent the voxels in 3D while the rows represent the different timepoints.
#' 
#' @title 4d nifti => 2d big matrix
#' @slot header A list that holds attributes of a nifti image
#' @slot mask A logical vector (represents the voxels in the original 4D image whose values are given in each column of \code{address}); its length must equal the product of the dim attribute in the header
#' @slot address Object of class \code{externalptr} points to the memory location of the C++ data structure.
setClass(
    Class = "big.nifti4d",
    representation = representation(),
    contains = c("header", "big.niftiXd")
)



