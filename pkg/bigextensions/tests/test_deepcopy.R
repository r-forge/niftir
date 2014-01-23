library(bigextensions)
library(testthat)
library(stringr)

context("Testing BeDeepCopy")

gen_data <- function(nr=100, nc=10, ...) {
    as.big.matrix(matrix(rnorm(nr*nc), nr, nc))
}

test_that("bedeepcopy gives the same results as deepcopy with no y", {
    x <- gen_data()
    
    ref.y <- deepcopy(x, rows=1:50, cols=1:5)
    comp.y <- bedeepcopy(x, x.rows=1:50, x.cols=1:5)
    
    expect_that(ref.y[,], equals(comp.y[,]))
})

test_that("bedeepcopy gives the same results as deepcopy with y", {
    x <- gen_data(shared=TRUE)
    y <- gen_data(shared=TRUE)
    
    ref.y <- deepcopy(y)
    sub.y <- sub.big.matrix(ref.y, firstRow=6, lastRow=55, firstCol=3, lastCol=7)
    sub.y <- deepcopy(x, rows=1:50, cols=1:5, y=sub.y)
    
    comp.y <- deepcopy(y)
    comp.y <- bedeepcopy(x, x.rows=1:50, x.cols=1:5, comp.y, y.rows=6:55, y.cols=3:7)
    
    expect_that(ref.y[,], equals(comp.y[,]))
})
