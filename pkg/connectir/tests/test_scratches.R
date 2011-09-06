library(connectir)
library(testthat)
library(stringr)

context("testing stuff")

# input for different tests
create_many_big_matrices <- function(nsubs=10, xdim=20, ydim=10, ...) {
    lapply(1:nsubs, function(x) as.big.matrix(scale(matrix(rnorm(xdim*ydim), xdim, ydim)), ...))
}

test_that("parts of compute_subdist_worker2 works", {
    nvoxs <- 100; nsubs <- 10
    funclist <- create_many_big_matrices(nsubs=nsubs, ydim=nvoxs)
    dmat <- big.matrix(nsubs^2, nvoxs, type="double", shared=FALSE)
    
    incols <- c(1,10)
    cormaps <- vbca_batch2(funclist, incols, ztransform=TRUE, shared=FALSE)
    
    i <- 2
    inds <- incols[1]:incols[2]
    voxs <- 1:nvoxs
    subsMap <- big.matrix(nvoxs-1, nsubs, type="double", shared=FALSE)
    
    # Combining and scaling
    .Call("CombineSubMapsMain", cormaps, subsMap@address, as.double(i), 
          as.double(voxs[-inds[i]]), as.double(nvoxs-1), as.double(nsubs))
    ref <- scale(sapply(cormaps, function(x) x[i,voxs[-inds[i]]]))
    expect_that(ref, is_equivalent_to(as.matrix(subsMap)))
    
    # Correlation
    big_cor(x=subsMap, z=dmat, z_firstCol=inds[i], z_lastCol=inds[i])
    m <- matrix(dmat[,inds[i]], nsubs, nsubs)
    mref <- cor(ref)
    expect_that(m, equals(mref))
    
    # Distance
    .Call("big_add_scalar", dmat, as.double(-1), as.double(1), 
            as.double(inds[i]), as.double(inds[i]));
    m <- matrix(dmat[,inds[i]], nsubs, nsubs)
    mref <- 1-cor(ref)
    expect_that(m, equals(mref))
})
