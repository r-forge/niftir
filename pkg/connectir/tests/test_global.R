library(connectir)
library(testthat)
library(stringr)
library(bigmemory)

context("global connectivity")

# input for different tests
create_big_matrix <- function(xdim=20, ydim=10, ...) {
    as.big.matrix(scale(matrix(rnorm(xdim*ydim), xdim, ydim)), ...)
}

# simpler function
simple_global <- function(bm) {
    nvoxs <- ncol(bm)
    cormat <- vbca(bm, 1:nvoxs)[,]
    gcor <- rowMeans(cormat)
    return(gcor)
}

# test rowsums
test_that("bm_rowmean works", {
    x <- create_big_matrix(xdim=200, ydim=100)
    ref <- rowMeans(x[,])
    comp <- bm_rowmean(x)
    expect_that(ref, equals(comp))
})

# test main function
test_that("gcor works", {
    x <- create_big_matrix(xdim=200, ydim=100)
    ref <- simple_global(x)
    comp <- gcor(x, blocksize=10)
    expect_that(ref, equals(comp))
})
