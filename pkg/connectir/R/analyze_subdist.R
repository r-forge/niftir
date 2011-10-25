
# Does SVM cross-validation for each voxel's subject distances matrix
# note that this won't save the fitted model, only returns the total accuracy of each cross-validation
svm_subdist_cross <- function(Xs, y, memlimit=2, bpath=NULL, 
                                cross=10, kernel="linear", 
                                verbose=T, parallel=F, ...) 
{
    vcat(verbose, "SVM on Subject Distances")
    
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
    
    # TODO: test this function
    fun1 <- function(dmat, y, cross, kernel, ...) {
        fit <- svm(dmat, y, cross=cross, kernel=kernel, fitted=FALSE, ...)
        fit$tot.accuracy
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    dmat <- matrix(0, nsubs, nsubs)
    for (bi in 1:blocks$n) {
        vcat(verbose, "...block %i", bi)
        first <- blocks$starts[bi]; last <- blocks$ends[bi]
        res <- llply(first:last, function(i) {
            X <- sub.big.matrix(Xs, firstCol=i, lastCol=i, backingpath=bpath)
            dcopy(N=nr, Y=dmat, X=X)
            fun1(dmat, y, cross, kernel, ...)
        }, .progress=progress, .inform=verbose, .parallel=parallel)
        accuracies <- c(accuracies, unlist(res))
        Xs <- free.memory(Xs, bpath)
        gc(FALSE, TRUE)
    }
    
    return(accuracies)
}

# TODO: same thing as above but for kmeans
