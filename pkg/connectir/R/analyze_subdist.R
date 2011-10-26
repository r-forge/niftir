base_analyze_subdist <- function(FUN, Xs, y, memlimit=2, bpath=NULL, 
                                verbose=T, parallel=F, ...) 
{
    vcat(verbose, "...checks/setup")
    if (!is.big.matrix(Xs))
        stop("input Xs must be big matrix")
    if (!is.shared(Xs))
        stop("input Xs must be a SHARED big matrix")
    if (is.filebacked(Xs) && is.null(bpath))
        stop("must specify bpath when Xs input is file-backed")
    if (!is.null(bpath) && !file.exists(bpath))
        vstop("couldn't find backing path '%s'", bpath)
    nvoxs <- ncol(Xs)
    nsubs <- length(y)
    nr <- nrow(Xs)
    if (nsubs != sqrt(nr))
        stop("length mismatch between 'y' labels and rows of Xs")
    
    vcat(verbose, "...setting memory limit of %f GB", memlimit)
    blocksize <- get_sdist_analysis_limit(memlimit, Xs)
    vcat(verbose, "...setting block size to %i (out of %i voxels)", blocksize, nvoxs)
    
    blocks <- niftir.split.indices(1, nvoxs, by=blocksize)
    progress <- ifelse(verbose, "text", "none")
    
    if (!is.null(bpath)) {
        vcat(verbose, "...freeing memory of subject distances")
        Xs <- free.memory(Xs, bpath)
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    dmat <- matrix(0, nsubs, nsubs)
    results <- c()
    for (bi in 1:blocks$n) {
        vcat(verbose, "...block %i", bi)
        first <- blocks$starts[bi]; last <- blocks$ends[bi]
        res <- llply(first:last, function(i) {
            X <- sub.big.matrix(Xs, firstCol=i, lastCol=i, backingpath=bpath)
            dcopy(N=nr, Y=dmat, X=X)
            FUN(dmat, y, ...)
        }, .progress=progress, .inform=verbose, .parallel=parallel)
        results <- c(results, unlist(res))
        Xs <- free.memory(Xs, bpath)
        gc(FALSE, TRUE)
    }
    
    return(results)
}

svm_subdist_cross <- function(Xs, y, memlimit=2, bpath=NULL, 
                                cross=10, kernel="linear", 
                                verbose=T, parallel=F, ...) 
{
    vcat(verbose, "SVM on Subject Distances")
        
    # TODO: test this function
    FUN <- function(dmat, y, cross, kernel, ...) {
        fit <- svm(dmat, y, cross=cross, kernel=kernel, fitted=FALSE, ...)
        fit$tot.accuracy
    }
    
    base_analyze_subdist(FUN, Xs, y, memlimit, bpath, verbose, parallel, 
                         cross=cross, kernel=kernel, ...)
}

kmeans_subdist_cross <- function(Xs, y, memlimit=2, bpath=NULL, 
                                    iter.max=200, nstart=20, 
                                    verbose=T, parallel=F, ...) 
{
    vcat(verbose, "K-means on Subject Distances")
    
    k <- length(unique(y))
    if (k != 2)
        stop("kmeans subdist can only be done right now with 2 groups")
    
    # TODO: test this function
    FUN <- function(dmat, y, k, iter.max, nstart, ...) {
        ks <- kmeans(dmat, k, iter.max, nstart, ...)$cluster
        acc <- sum(diag(table(ks, y)))/length(y)
        if (acc < 0.5)  acc <- 1 - acc
        return(acc)
    }
    
    base_analyze_subdist(FUN, Xs, y, memlimit, bpath, verbose, parallel, 
                         k=k, iter.max=iter.max, nstart=nstart, ...)
}


# TODO: same thing as above but for kmeans
