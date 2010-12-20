\name{as.nifti4d}
\alias{as.nifti4d}
\alias{as.nifti4d}
\title{Create a nifti4d object...}
\usage{as.nifti4d(x)   # when x is a nifti object or big.nifti4d object
as.nifti4d(x, header)   # when x is an array
as.nifti4d(x, header, mask) # when x is a matrix}
\description{Create a nifti4d object}
\author{Zarrar Shehzad}
\value{\code{nifti4d} object}
\seealso{\code{\link{as.nifti}}, \code{\link{as.nifti4d-methods}}}
\keyword{methods}
\arguments{\item{x}{4d \code{nifti}, 4d \code{array}, \code{big.nifti4d}, or \code{matrix}}
\item{header}{list of header attributes (required when \code{x} is an
\code{array} or \code{matrix})}
\item{mask}{logical vector specifying which voxels from 3D image are
specified with \code{x} (only required when \code{x} is a \code{matrix})}
}
\examples{as.nifti4d(array(0, c(10,10,10,10)), create.header())
as.nifti4d(nifti(0, dim=c(10,10,10,10)))   # should give same thing as above}

