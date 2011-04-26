library('testthat')

test.vbca <- function() {
    require(connectir)
    
    nvoxs <- 1000
    ntpts <- 197
    
    mat <- matrix(rnorm(nvoxs*ntpts), ntpts, nvoxs)
    bm <- as.big.matrix(mat)
    scale(bm, to.copy=FALSE)
    
    tmp1 <- vbca(bm, 10:15)
    tmp2 <- cor(mat[,1:15])
    
    expect_equal(as.vector(tmp1[,10:15]), as.vector(tmp2[10:15,10:15]))
}

test.batch.vbca <- function() {
    require(connectir)

    nvoxs <- 1000
    nsubs <- 5
    ntpts <- 197
    nseeds <- 10
    
    subs.mats <- lapply(1:nsubs, function(i) matrix(rnorm(nvoxs*ntpts), ntpts, nvoxs))
    subs.bigmats <- lapply(subs.mats, function(x) as.big.matrix(scale(x)))
    subdist <- big.matrix(nsubs^2, nseeds)
    
    system.time(cormaps_list <- vbca_batch(subs.bigmats, 4:8))
    system.time(comp <- lapply(1:nsubs, function(i) cor(subs.bigmats[[i]][,4:8], subs.bigmats[[i]][,])))
    
    expect_equal(as.vector(cormaps_list[[1]][,]), as.vector(comp[[1]]))
}


test.combine.submaps <- function() {
    require(connectir)

    nvoxs <- 100
    nsubs <- 5
    nseeds <- 10
    
    subs.mats <- lapply(1:nsubs, function(i) matrix(rnorm(nvoxs*nseeds), nseeds, nvoxs))
    subs.bigmats <- lapply(subs.mats, function(x) as.big.matrix(x))

    out.bigmat <- big.matrix(nvoxs, nsubs)    
    .Call("CombineSubMapsMain", subs.bigmats, out.bigmat@address, as.double(1), as.double(1:nvoxs), as.double(nvoxs), as.double(nsubs))
    out.mat <- sapply(1:nsubs, function(i) subs.mats[[i]][1,])
    out.mat <- scale(out.mat)
    expect_equal(as.vector(out.mat), as.vector(out.bigmat[,]))
    
    out.bigmat <- big.matrix(nvoxs, nsubs)    
    .Call("CombineSubMapsMain", subs.bigmats, out.bigmat@address, as.double(3), as.double(1:nvoxs), as.double(nvoxs), as.double(nsubs))
    out.mat <- sapply(1:nsubs, function(i) subs.mats[[i]][3,])
    out.mat <- scale(out.mat)  
    expect_equal(as.vector(out.mat), as.vector(out.bigmat[,]))   
}

compare.subdist <- function(sub.cormaps, cols, outmat=NULL, FUN=cor, ...) {
    sapply(1:length(cols), function(i) {
        as.vector(1 - FUN(sapply(sub.cormaps, function(x) x[i,-cols[i]]), ...))
    })
}

test.compute.subdist <- function() {
    require(connectir)

    nvoxs <- 1000
    nsubs <- 5
    ntpts <- 197
    nseeds <- 10
    
    subs.mats <- lapply(1:nsubs, function(i) matrix(rnorm(nvoxs*ntpts), ntpts, nvoxs))
    subs.bigmats <- lapply(subs.mats, function(x) as.big.matrix(scale(x)))
    
    subdist <- big.matrix(nsubs^2, nseeds)
    compute_subdist(subs.bigmats, subdist, seed_inds=1:nseeds, blocksize=4, ztransform=FALSE, start=1, verbose=TRUE)
    
    cormaps_compare <- lapply(1:nsubs, function(i) cor(subs.bigmats[[i]][,1:nseeds], subs.bigmats[[i]][,]))
    subdist_compare <- compare.subdist(cormaps_compare, 1:nseeds)
    
    expect_equal(as.vector(subdist[,]), as.vector(subdist_compare))
    
    
    subdist <- big.matrix(nsubs^2, nseeds)
    compute_subdist(subs.bigmats, subdist, seed_inds=5:nseeds, blocksize=4, ztransform=FALSE, start=1, verbose=TRUE)
    
    cormaps_compare <- lapply(1:nsubs, function(i) cor(subs.bigmats[[i]][,5:nseeds], subs.bigmats[[i]][,]))
    subdist_compare <- compare.subdist(cormaps_compare, 5:nseeds)
    
    expect_equal(as.vector(subdist[,1:6]), as.vector(subdist_compare))
}


test.compute.subdist2 <- function() {
    require(connectir)
    
    nvoxs <- 1000
    nsubs <- 5
    ntpts <- 197
    nseeds <- 10
    
    subs.mats <- lapply(1:nsubs, function(i) matrix(rnorm(nvoxs*ntpts), ntpts, nvoxs))
    funclist <- lapply(subs.mats, function(x) as.big.matrix(scale(x)))
    subdist <- big.matrix(nsubs^2, nseeds)
    
    seed_inds <- 1:10
    ztransform <- FALSE
    verbose <- TRUE
    start <- 1
    nseeds <- length(seed_inds)
    blocksize <- 4
    
    
    nseeds <- length(seed_inds)
    blocks <- niftir.split.indices(start, nseeds, by=blocksize)
    
    dfun <- function(i, blocks, seed_inds, funclist, subdist, ztransform, verbose, pb) {
        if (verbose)
            update(pb, i)
        inds_CHUNK <- seed_inds[blocks$starts[i]:blocks$ends[i]]
        cormaps_list <- vbca_batch(funclist, inds_CHUNK, ztransform=ztransform)
        subdist_CHUNK <- sub.big.matrix(subdist, firstCol=blocks$starts[i], lastCol=blocks$ends[i])
        compute_subdist_worker(cormaps_list, inds_CHUNK, subdist_CHUNK)
    }
    
    if (verbose)
        pb <- progressbar(blocks$n)
    else
        pb <- NULL
    
    if (getDoParRegistered()) {
        foreach(i=1:blocks$n, .packages=c("connectir")) %dopar% 
            dfun(i, blocks, seed_inds, funclist, subdist, ztransform, verbose, pb)
    }
    else {
        foreach(i=1:blocks$n, .packages=c("connectir")) %do% 
            dfun(i, blocks, seed_inds, funclist, subdist, ztransform, verbose, pb)
    }
    
    if (verbose)
        end(pb)
        
    tmpinds <- blocks$starts[1]:blocks$ends[1]
    cormaps_compare <- lapply(1:nsubs, function(i) cor(funclist[[i]][,tmpinds], funclist[[i]][,]))
    subdist_compare <- compare.subdist(cormaps_compare, tmpinds)
    
    all.equal(as.vector(cormaps_compare[[2]][,]), as.vector(cormaps_list[[2]][,])) 
    
}


timing.combine.submaps <- function() {
    require(connectir)

    nvoxs <- 10000
    nsubs <- 5
    nseeds <- 10
    
    subs.mats <- lapply(1:nsubs, function(i) matrix(rnorm(nvoxs*nseeds), nseeds, nvoxs))
    subs.bigmats <- lapply(subs.mats, function(x) as.big.matrix(x))
    subdist <- big.matrix(nsubs^2, nseeds)
    
    system.time(tmp1 <- compare.subdist(subs.bigmats, 1:nseeds))
    system.time(tmp2 <- compute_subdist_worker(subs.bigmats, 1:nseeds))
    system.time(compute_subdist_worker(subs.bigmats, 1:nseeds, subdist))
    
    expect_equal(as.vector(tmp1[,]), as.vector(tmp2[,]))
    
    print(system.time(out.bigmat <- big.matrix(nvoxs, nsubs)))
    print(system.time({
        .Call("CombineSubMapsMain", subs.bigmats, out.bigmat@address, as.double(1), as.double(nvoxs), as.double(nsubs))
    }))
    print(system.time({
        out.mat <- sapply(1:nsubs, function(i) subs.mats[[i]][,1])
        out.mat <- scale(out.mat)
    }))
    
}

test.mdmr.prepare.model <- function() {
    model <- data.frame(group=factor(c(rep(c(1,2,3),30))), age=rnorm(90), random=rnorm(90))
    return(mdmr.prepare.model(.~group+age+random, model))
}

test.mdmr <- function() {
    nsubs <- 90
    nvoxs <- 120
    x <- as.big.matrix(matrix(rnorm(nsubs^2 * nvoxs), nsubs^2, nvoxs))
    model <- data.frame(group=factor(c(rep(1:3,nsubs/3))), age=rnorm(nsubs), random=rnorm(nsubs))
    return(mdmr(x, .~group+age+random, model, nperms=99, factors.to.perm=c("group"), block.size=50))
}

test.subdist <- 