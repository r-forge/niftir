vox_base_cross <- function(FUN, funclist1, y, blocksize, 
                              funclist2=funclist1, 
                              verbose=TRUE, parallel=FALSE, shared=parallel, 
                              ztransform=FALSE, recursive.unlist = TRUE, 
                              ...)
{
    vcat(verbose, "...checks")
    if (length(funclist1) != length(funclist2))
        stop("length mismatch between first and second set of functionals")
    if (length(funclist1) != length(y))
        stop("length mismatch between first set of functionals and labels")
    
    vcat(verbose, "...setup")
    nsubs <- length(funclist1)
    nvoxs1 <- ncol(funclist1[[1]])
    nvoxs2 <- ncol(funclist2[[1]])
    voxs1 <- 1:nvoxs1
    voxs2 <- 1:nvoxs2
    blocks <- niftir.split.indices(1, nvoxs1, by=blocksize)
    progress <- ifelse(verbose, "text", "none")
    
    gfun <- function(bi) {
        vcat(verbose, "...block %i", bi)
        
        first <- blocks$starts[bi]; last <- blocks$ends[bi]
        sub_voxs <- first:last
        sub_nvoxs <- length(sub_voxs)
                
        # correlate
        vcat(verbose, "....correlate")
        subs.cormaps <- vbca_batch3(funclist1, c(first, last), funclist2, 
                                    ztransform=ztransform, 
                                    type="double", shared=shared)
        
        # loop through each seed
        vcat(verbose, '....classify')
        sfun <- function(si) {
            # combine and transpose correlation maps for 1 seed
            seedCorMaps <- big.matrix(nsubs, nvoxs2, type="double", shared=FALSE)
            .Call("subdist_combine_and_trans_submaps", subs.cormaps, as.double(si), 
                  as.double(voxs2), seedCorMaps, PACKAGE="connectir")
            #X <- matrix(NA, nsubs, nvoxs2)
            #dcopy(X=seedCorMaps, Y=X)
            X <- seedCorMaps[,]
            
            # some classification?
            res <- FUN(X, y, ...)
            
            rm(seedCorMaps, cmaps); gc(FALSE, TRUE)
            
            return(res)
        }
        
        res <- llply(1:sub_nvoxs, sfun, .parallel=parallel, .progress=progress, .inform=T)
                
        rm(subs.cormaps); gc(FALSE, TRUE)
        
        return(res)
    }
    
    vcat(verbose, "...%i blocks", blocks$n)
    res <- lapply(1:blocks$n, gfun)
    res <- unlist(res, recursive=recursive.unlist)
        
    return(res)
}

vox_svm_cross <- function(funclist1, y, blocksize, 
                          funclist2=funclist1, 
                          cross=10, kernel="linear", 
                          verbose=T, parallel=F, shared=parallel, 
                          ztransform=TRUE, 
                          ...) 
{
    vcat(verbose, "SVM on Subject Distances")
    library(e1071)
        
    # TODO: test this function
    FUN <- function(X, y, cross, kernel, ...) {
        fit <- svm(X, y, cross=cross, kernel=kernel, fitted=FALSE, ...)
        res <- c(fit$tot.accuracy, fit$tot.MSE[[1]])
        res[!is.null(res)]
    }
    
    tmp <- vox_base_cross(FUN, funclist1, y, blocksize, funclist2, verbose, 
                            parallel, shared, ztransform, recursive.unlist=TRUE,  
                            cross=cross, kernel=kernel, ...)
    
    return(tmp)
}
