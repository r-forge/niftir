library(testthat)
library(stringr)
library(bigmemory)

context("kendall")

# input for different tests
create_many_big_matrices <- function(nsubs=10, xdim=20, ydim=10, ...) {
    lapply(1:nsubs, function(x) as.big.matrix(scale(matrix(rnorm(xdim*ydim), xdim, ydim)), ...))
}

# simpler function
simple_kendall <- function(bigmats) {
    nsubs <- length(bigmats)
    nvoxs <- ncol(bigmats[[1]])
    cormats <- vbca_batch(bigmats, 1:nvoxs)
    coeffs <- sapply(1:nvoxs, function(i) {
        xmat <- sapply(cormats, function(x) x[,i])
        kendall_ref(xmat)
    })
    return(coeffs)
}

# test
test_that("works", {
    xs <- create_many_big_matrices(nsub=10, xdim=200, ydim=100)
    ref <- simple_kendall(xs)
    comp <- kendall(xs, blocksize=10)
    expect_that(ref, equals(comp))
})
