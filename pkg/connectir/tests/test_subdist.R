library(connectir)
library(testthat)
library(stringr)

context("Subject Distances")

create_data <- function(n=10, ncols=5) {
    xm <- matrix(runif((n^2)*ncols, min=0, max=2), n^2, ncols)
    xbm <- as.big.matrix(xm, type="double", shared=FALSE)
    return(list(m=xm, bm=xbm))
}

# G <- -0.5 * dmat^2 %*% (I - ones %*% t(ones)/n)
test_that("creation of gower matrix works", {
    dat <- create_data()
    
    n <- sqrt(nrow(dat$m))
    I <- diag(n)
    ones <- matrix(1, n, 1)
    adj <- I - tcrossprod(ones)/n
    ref <- apply(dat$m, 2, function(x) {
        dmat <- matrix(x, n, n)
        A <- -1*(dmat^2)/2
        as.vector(A %*% adj)
    })
    
    G <- big.matrix(nrow(dat$bm), ncol(dat$bm), type="double", shared=FALSE)
    nc <- ncol(dat$bm)
    comp <- .Call("big_gower", dat$bm, G, as.double(1), as.double(nc), 
                  as.double(1), as.double(nc), PACKAGE="connectir")
    expect_that(ref, is_equivalent_to(as.matrix(G)))
})
