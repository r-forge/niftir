\name{as.nifti3d}
\alias{as.nifti3d}
\alias{as.nifti3d}
\title{Create a nifti3d object...}
\usage{as.nifti3d(x)   # when x is a nifti object
as.nifti3d(x, header)   # when x is an array
as.nifti3d(x, header, mask) #when x is a vector}
\description{Create a nifti3d object}
\author{Zarrar Shehzad}
\value{\code{nifti4d} object}
\seealso{\code{\link{as.nifti}}, \code{\link{as.nifti3d-methods}}}
\keyword{methods}
\arguments{\item{x}{3D \code{nifti}, 3D \code{array}, or \code{vector}}
\item{header}{list of header attributes (required when \code{x} is an
\code{array} or \code{vector})}
\item{mask}{logical vector specifying which voxels from 3D image are
specified with \code{x} (only required when \code{x} is a \code{vector})}
}
\examples{as.nifti3d(array(0, c(10,10,10)), create.header())
as.nifti3d(nifti(0, dim=c(10,10,10)))   # should give same thing as above}

