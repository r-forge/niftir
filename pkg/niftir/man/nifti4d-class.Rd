\name{nifti4d-class}
\alias{nifti4d-class}
\alias{nifti4d}
\title{4d nifti => 2d matrix}
\description{Represents a 4D nifti image as a 2D matrix}
\details{Assuming that the 4D image is akin to 3D+time, then the columns represent the voxels in 3D while the rows represent the different timepoints.}
\section{Slots}{\describe{\item{\code{header}:}{}\item{\code{mask}:}{}\item{\code{.Data}:}{}}}
\section{Extends}{\code{\linkS4class{header}}, \code{\linkS4class{niftiXd}}, \code{\linkS4class{matrix}}}

