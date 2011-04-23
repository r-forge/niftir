#' @nord
setGeneric("scale",
    function(x, center=TRUE, scale=TRUE, ...)
        standardGeneric('scale')
)

#' @nord
setMethod("scale",
    signature(x="big.matrix"),
    function(x, center=TRUE, scale=TRUE, to.copy=TRUE) {
        # Copy
        if (to.copy)
            y <- zarrar.deepcopy(x)
        else
            y <- x

        # Get mean and/or sd column values
        # and center and/or scale
        if (center) {
            xmean <- colmean(y, na.rm=TRUE)
            .Call("BigSubtractColsMain", y@address, as.double(xmean), TRUE)
        }
        if (scale) {
            xsd <- colsd(y, na.rm=TRUE)
            .Call("BigDivideColsMain", y@address, as.double(xsd), TRUE)
        }
        
        return(y)
    }
)

#' @nord
setMethod("scale",
    definition=function(x, center=TRUE, scale=TRUE)
        scale.default(x, center, scale)
)

#' @nord
setGeneric("do.operator",
    function(x, y, operator, z, ...)
        standardGeneric('do.operator')
)

#' @nord
setMethod("do.operator",
    signature(x="big.matrix", y="big.matrix", operator="character", z="big.matrix"),
    function(x, y, operator, z, xCols=1:ncol(x), yCols=1:ncol(y), zCols=1:ncol(z), xRows=1:nrow(x), yRows=1:nrow(y), zRows=1:nrow(z)) {
        if (!all(length(xCols)==length(yCols)) && !all(length(xCols)==length(zCols)))
            stop("Length of column indices must be the same across matrices")
        if (!all(length(xRows)==length(yRows)) && !all(length(xRows)==length(zRows)))
            stop("Length of row indices must be the same across matrices")
        switch(operator,
            '+' = .Call("BigAddMatsMain", x@address, y@address, z@address, as.double(xCols), as.double(yCols), as.double(zCols), as.double(xRows), as.double(yRows), as.double(zRows)),
            '-' = .Call("BigSubtractMatsMain", x@address, y@address, z@address, as.double(xCols), as.double(yCols), as.double(zCols), as.double(xRows), as.double(yRows), as.double(zRows)),
            '*' = .Call("BigMultiplyMatsMain", x@address, y@address, z@address, as.double(xCols), as.double(yCols), as.double(zCols), as.double(xRows), as.double(yRows), as.double(zRows)),
            '/' = .Call("BigDivideMatsMain", x@address, y@address, z@address, as.double(xCols), as.double(yCols), as.double(zCols), as.double(xRows), as.double(yRows), as.double(zRows)),
            stop("operator must be +, -, *, or /")
        )
        return(NULL)
    }
)

#' @nord
# NOTE: THIS IS CURRENTLY SPECIFIC TO NEUROCOR/6D OBJECTS...BAD
setMethod("do.operator",
    signature(x="big.matrix", y="big.matrix", operator="character", z="missing"),
    function(x, y, operator, ...) {
        zlist <- niftir.big.matrix(nrow(x), ncol(x), ...)
        do.operator(x, y, operator, zlist$bm)
        if (is.big.niftiXd(x)) {
            # create new class
            args <- list(class(x))
            zslots <- lapply(slotNames(x), function(sn) slot(x, sn))
            names(zslots) <- slotNames(x)
            zslots$address <- zlist$bm@address
            zslots$backingfile <- zlist$bf
            zslots$descriptorfile <- zlist$df
            args <- append(args, zslots)
            z <- do.call("new", args)
            describe.ncor <- list(descriptorfile=args$descriptorfile, backingfile=args$backingfile, header=args@header, mask=args@mask, rowmask=args@mask)
            fname <- sprintf("%s.neurocor", rmext(args$backingfile))
            dput(describe.ncor, file=fname)
        } else {
            z <- zlist$bm
        }
        return(z)
    }
)

# Calculate the sum across many big matrices
setGeneric("sum.bigmats",
    function(x, blocksize=1000, verbose=TRUE, ...)
        standardGeneric('sum.bigmats')
)

setMethod("sum.bigmats",
    signature(x="list"),
    function(x, blocksize=1000, colInds=NULL, rowInds=NULL, verbose=T, toSquare=F, backingpaths=NULL) {
        
        if (is.big.nifti4d(x[[1]]))
            backingpaths <- sapply(x, function(xx) dirname(xx@backingfile))
        if (!is.null(backingpaths) && length(x) != length(backingpaths))
            stop("list length between x and backingpaths input must be the same")
        
        nx <- length(x)
        # Check that lengths match
        if (!is.null(colInds) && length(x) != length(colInds))
            stop("list length between x and colInds input must be the same")
        if (!is.null(rowInds) && length(x) != length(rowInds))
            stop("list length between x and rowInds input must be the same")
        # Check that all objects the same
        out.type <- typeof(x[[1]])
        if (!all(out.type == sapply(x[-1], typeof)))
            stop("All input objects must have the same data type")
        # Check columns
        if (is.null(colInds)) {
            refCol <- ncol(x[[1]])
            colInds <- lapply(1:nx, function(i) {
                iCol <- ncol(x[[i]])
                if (iCol != refCol)
                    stop("All inputs must have the same column dimensions")
                iCol
            })
        } else {
            refCol <- length(colInds[[1]])
            lapply(1:nx, function(i) {
                if (length(colInds[[1]]) != refCol)
                    stop("All inputs must have the same number of column indices")
            })
        }
        # Check rows
        if (is.null(rowInds)) {
            refRow <- nrow(x[[1]])
            rowInds <- lapply(1:nx, function(i) {
                iRow <- nrow(x[[i]])
                if (iRow != refRow)
                    stop("All inputs must have the same row dimensions")
                iRow
            })
        } else {
            refRow <- length(rowInds[[1]])
            lapply(1:nx, function(i) {
                if (length(rowInds[[1]]) != refRow)
                    stop("All inputs must have the same number of row indices")
            })
        }
        
        # Create a new bigmat object
        outmat <- big.matrix(refRow, refCol, type=out.type, init=0)
        
        # Sum Function
        if (toSquare) {
            # Sum Squared Function
            sfun <- function(j) {
                # Copy over to memory
                tmp.mat <- zarrar.deepcopy(x[[xi]], cols=colInds[[xi]][j], rows=rowInds[[xi]])

                # Square
                .Call("BigPowMain", tmp.mat@address, PACKAGE="niftir")

                # Add Together
                do.operator(tmp.mat, outmat, "+", outmat, yCols=j, zCols=j)
            }
        } else {
            sfun <- function(j) {
                do.operator(x[[xi]], outmat, "+", outmat, xCols=colInds[[xi]][j], xRows=rowInds[[xi]], yCols=j,  zCols=j)
            }
        }
        
        # Loop through elements adding to output
        blocks <- niftir.split.indices(1, refCol, by=blocksize)
        #isFB <- sapply(x, function(xx) is.filebacked(xx))
        #anyFB <- any(isFB)
        if (verbose)
            pb <- progressbar(blocks$n)
        for (bi in 1:blocks$n) {
            if (verbose) {
                update(pb, bi)
                stime <- proc.time()[3]
            }
            
            for (xi in 1:nx) {
                # Break up block into sub-blocks
                if (getOption("neurocor.allow.parallel") && getDoParRegistered())
                    foreach(j=blocks$starts[bi]:blocks$ends[bi]) %dopar% sfun(j)
                else
                    sfun(blocks$starts[bi]:blocks$ends[bi])
                x[[xi]] <- free.memory(x[[xi]], backingpath=backingpaths[[xi]])
            }
            
            if (verbose) {
                etime <- proc.time()[3]
                totaltime <- ((etime-stime)*(blocks$n-bi))/60 # minutes
                cat("\nEstimated completion time:", totaltime, "mins\n")
            }
        }
        
        return(outmat)
    }
)

# Calculate the mean across many big matrices
setGeneric("mean.bigmats",
    function(x, ...)
        standardGeneric('mean.bigmats')
)

setMethod("mean.bigmats",
    signature(x="list"),
    function(x, ...) {
        outmat <- sum.bigmats(x, ...)
        .Call("BigDivideScalarMain", outmat@address, as.double(length(x)), TRUE, PACKAGE="niftir")
        return(outmat)
    }
)

# Calculate variance across many big matrices
setGeneric("var.bigmats",
    function(x, mx, ...)
        standardGeneric('var.bigmats')
)

setMethod("var.bigmats",
    signature(x="list", mx="big.matrix"),
    function(x, mx, blocksize=1000, colInds=NULL, rowInds=NULL, verbose=T, mx.backingpath=NULL) {
        
        # Check that lengths match
        if (!is.null(colInds) && length(x) != length(colInds))
            stop("list length between x and colInds input must be the same")
        if (!is.null(rowInds) && length(x) != length(rowInds))
            stop("list length between x and rowInds input must be the same")
        
        nx <- length(x)
        # Check that all objects the same
        out.type <- typeof(x[[1]])
        if (!all(out.type == sapply(x[-1], typeof)))
            stop("All input objects must have the same data type")
        # Check columns
        if (is.null(colInds)) {
            refCol <- ncol(mx)
            colInds <- lapply(1:nx, function(i) {
                iCol <- ncol(x[[i]])
                if (iCol != refCol)
                    stop("All inputs must have the same column dimensions")
                iCol
            })
        } else {
            refCol <- ncol(mx)
            lapply(1:nx, function(i) {
                if (length(colInds[[1]]) != refCol)
                    stop("All inputs must have the same number of column indices")
            })
        }
        # Check rows
        if (is.null(rowInds)) {
            refRow <- nrow(mx)
            rowInds <- lapply(1:nx, function(i) {
                iRow <- nrow(x[[i]])
                if (iRow != refRow)
                    stop("All inputs must have the same row dimensions")
                iRow
            })
        } else {
            refRow <- nrow(mx)
            lapply(1:nx, function(i) {
                if (length(rowInds[[1]]) != refRow)
                    stop("All inputs must have the same number of row indices")
            })
        }
        
        # Create a new bigmat object
        outmat <- big.matrix(refRow, refCol, type=out.type, init=0)
        
        # Var Function
        sfun <- function(j) {
            # Create Temporary Big Matrix
            tmp.mat <- big.matrix(refRow, length(j), type=out.type, init=0)
            
            # De-mean [x-mean(x)]
            do.operator(x[[xi]], mx, "-", tmp.mat, xCols=colInds[[xi]][j], xRows=rowInds[[xi]], yCols=j)
            
            # Square [(x-mean(x))^2]
            .Call("BigPowMain", tmp.mat@address, PACKAGE="niftir")
            
            # Add Together
            do.operator(tmp.mat, outmat, "+", outmat, yCols=j, zCols=j)
            
            rm(tmp.mat)
            gc()
        }
        
        # Loop through elements adding to output
        blocks <- niftir.split.indices(1, refCol, by=blocksize)
        if (verbose)
            pb <- progressbar(blocks$n)
        for (bi in 1:blocks$n) {
            if (verbose) {
                update(pb, bi)
                stime <- proc.time()[3]
            }
            
            for (xi in 1:nx) {
                # Break up block into sub-blocks
                if (getOption("neurocor.allow.parallel") && getDoParRegistered())
                    foreach(j=blocks$starts[bi]:blocks$ends[bi]) %dopar% sfun(j)
                else
                    sfun(blocks$starts[bi]:blocks$ends[bi])
                x[[xi]] <- free.memory(x[[xi]])
            }
            mx <- free.memory(mx, mx.backingpath)
                        
            if (verbose) {
                etime <- proc.time()[3]
                totaltime <- ((etime-stime)*(blocks$n-bi))/60 # minutes
                cat("\nEstimated completion time:", totaltime, "mins\n")
            }
        }
        
        # Divide by df
        if (verbose)
            cat("Almost done...\n")
        .Call("BigDivideScalarMain", outmat@address, as.double(length(x)-1), TRUE, PACKAGE="niftir")
        
        return(outmat)
    }
)

