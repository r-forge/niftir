\name{read.nifti3d}
\alias{read.nifti3d}
\title{Read in a nifti3d object from a file...}
\usage{read.nifti3d(fname)}
\description{Read in a nifti3d object from a file}
\author{Zarrar Shehzad}
\value{\code{nifti4d} object}
\seealso{\code{\link{read.nifti}}, \code{\link{as.nifti3d}}}
\keyword{methods}
\arguments{\item{fname}{character specifying path to analyze/nifti file}
}
\examples{fname <- file.path(system.file("data", package="niftir"), "test.nii.gz") # 2x2x2 size file
# or fname <- file.choose()
nim <- read.nifti3d(fname)}

