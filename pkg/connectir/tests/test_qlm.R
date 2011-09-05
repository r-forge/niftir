library(connectir)
library(testthat)
library(stringr)

context("qlm")

create_data <- function(n=100, k=2, m=3) {
    y <- matrix(rnorm(n*m), n, m)
    X <- cbind(rep(1, n), matrix(rnorm(n*k), n, k))
    return(list(y=y, X=X))
}

test_that("qlm dd works?", {
    dat <- create_data()
    
    ref <- diag(solve(t(dat$X) %*% dat$X))      # will return a vector
    comp <- qlm_dd(as.big.matrix(dat$X))                # will return 1 column
    
    expect_that(ref, is_equivalent_to(comp))
})

test_that("qlm fit works?", {
    dat <- create_data()
    
    ref <- lm(dat$y ~ dat$X - 1)
    comp <- qlm_fit(as.big.matrix(dat$y), as.big.matrix(dat$X))
    
    expect_that(as.vector(ref$coef), equals(as.vector(comp$coef[,])))
    
    expect_that(as.vector(ref$residuals), equals(as.vector(comp$residuals[,])))
    
    mse <- colSums(ref$residuals^2)/ref$df.residual
    expect_that(mse, is_equivalent_to(comp$mse[,]))
})

test_that("qlm residuals works?", {
    dat <- create_data(k=2)
    
    fit <- lm(dat$y ~ dat$X - 1)
    residuals <- qlm_residuals(as.big.matrix(dat$y), as.big.matrix(dat$X))
    
    ref <- as.matrix(fit$residuals)
    comp <- as.matrix(residuals)
    expect_that(ref, is_equivalent_to(comp))
})

test_that("qlm contrasts works?", {
    dat <- create_data(k=2)
    
    fit <- qlm_fit(as.big.matrix(dat$y), as.big.matrix(dat$X))

    dd <- qlm_dd(as.big.matrix(dat$X))
    cons <- matrix(c(0,1,0,0,0,1), 2, 3, byrow=T)
    comp <- qlm_contrasts(fit, cons, dd)
    
    ref <- list()
    ref$coefficients <- cons %*% as.matrix(fit$coef)
    c_dd <- cons %*% dd
    ref$standard_errors <- sqrt(sapply(as.vector(fit$mse[,]), function(x) c_dd * x))
    ref$tvals <- ref$coefficients/ref$standard_errors
    
    expect_that(ref$coef, is_equivalent_to(as.matrix(comp$coef)))
    expect_that(ref$st, is_equivalent_to(as.matrix(comp$st)))
    expect_that(ref$tval, is_equivalent_to(as.matrix(comp$tval)))
})

