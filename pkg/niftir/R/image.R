# TODO: add image function (might want plot function for nifti4D)
# provide an underlay and an overlay
# converting slice from overlay to underlay
# t(y@header$qto.xyz) %*% t(x@header$qto.ijk)
# need make sure that off diagonal of 1:3,1:3 of above is 0 for the slice selector to be accurate
# c(10,1) %*% y@header$qto.xyz[1,c(1,4)] ... here can get the conversion for one slice

#' @nord
"+.nimage" <- function(p1, p2) {
    if (is.nifti(p2)) {
        p1 <- p1 + ni_overlay(p2)
        # p1$overlays[[length(p1$overlays)+1]] <- p2 # CHANGE THIS FORMAT
    } else {
        p <- switch(class(p2),
            list = {
                for (p3 in p2)
                    p1 <- p1 + p3
                p1
            },
            ni_slice = {
                p1$slice <- defaults(p2, p1$slice)
                p1
            },
            ni_layout = {
                p1$layout <- defaults(p2, p1$layout)
                p1
            },
            ni_underlay = {
                p1$underlay <- p2$underlay
                p1$underlay.opts <- p2$underlay.opts
                p1
            },
            ni_overlay = {
                p1$overlays <- append(p1$overlays, p2$overlays)
                p1$overlays.opts <- append(p1$overlays.opts, p2$overlays.opts)
                p1
            },
            stop("Class ", class(p2), " not recognized in +.nimagine")
        )
    }
        
    return(p1)
}

#' @nord
"%+%" <- `+.nimage`

#' @nord
nimage <- function(data, ...) UseMethod("nimage")

#' @nord
nimage.default <- function(...) {
    nimage(nifti(0), is.underlay.empty=T, ...)
}

#' @nord
nimage.list <- function(data, ...) {
    nimage(nifti(0), overlays=data, is.underlay.empty=T, ...)
}

#' @nord
nimage.nifti <- function(data, overlays=NULL, is.underlay.empty=F, ...) {
    if (is.underlay.empty)
        underlay <- NULL
    else
        underlay <- data
    
    p <- structure(
        list(
            underlay = NULL,
            underlay.opts = NULL,
            overlays = list(),
            overlays.opts = list()
        ),
        class="nimage"
    )
        
    if (!is.null(underlay))
        p <- p + ni_underlay(underlay, ...)
    
    if (!is.null(overlays))
        p <- p + ni_overlay(overlays)
    
    return(p)
}

# slice stuff
#' @nord
ni_slice <- function(slices=-0.5, orientation="sagittal", nslices=NULL, show.numbers=F) {
    if (!is.null(nslices))
        slices <- as.character(nslices)
    structure(list(
        slices=slices,
        orientation=orientation,
        show.numbers=show.numbers
    ), class="ni_slice")
}

# plot layout
#' @nord
ni_layout <- function(nc=NULL) {
    structure(list(
        nc = nc
    ), class="ni_layout")
}

# underlay
#' @nord
ni_underlay <- function(underlay=NULL, use.standard=NULL, col=c("white", colorRampPalette(c("black", "grey", "white"))(256)), ...) {
    # Use a standard brain/head
    ## get option
    use.standard = match(use.standard, c(NULL, "1b", "1h", "2b", "2h"))
    if (is.na(use.standard))
        stop("options for use.standard are 1b, 1h, 2b, and 2h. The numbers correspond to the resolution of the MNI152 standard while the letters correspond to whether it is a brain or head")
    ## get standard
    if (!is.null(use.standard)) {
        if (!is.null(underlay))
            warning("specified both an underlay and a use.standard option...going with the use.standard (eek)")
        basepath <- system.file("data", package="niftir")
        fname <- switch(use.standard,
            "1b" = "standard_brain_1mm.nii.gz",
            "1h" = "standard_head_1mm.nii.gz",
            "2b" = "standard_brain_2mm.nii.gz",
            "2h" = "standard_head_2mm.nii.gz"
        )
        underlay <- read.nifti(file.path(basepath, fname))
    }
    
    if (!is.null(underlay) && is.character(underlay))
        underlay <- read.nifti(underlay)
    
    if (!is.null(underlay) && !is.nifti(underlay))
        underlay <- as.nifti(underlay, ...)
    
    structure(list(
        underlay=underlay,
        underlay.opts=list(
            col=col # maybe will add more options later (for now this is a loner)
        )
    ), class="ni_underlay")
}

# overlays
## ex: thresh="overlay<2.3 | overlay>-2.3"
#' @nord
ni_overlay <- function(overlays, col=colorRampPalette(c("yellow", "orange", "red"))(256), colrange=NULL, thresh=NULL, to.binarize=FALSE, auto.colrange=TRUE, auto.thresh=TRUE, show.all.values=TRUE, ...) {
    tmp <- vector("list", length(overlays))
    outs <- structure(tmp, class="ni_overlay")
    overlays.opts <- list(col=col)

    # Save the overlays
    if (is.list(overlays)) {
        for (o in 1:length(overlays)) {
            overlay <- overlays[[o]]
            if (is.character(overlay))
                overlay <- read.nifti(overlay)
            outs[[o]]$overlays.opts <- overlays.opts
            outs[[o]]$overlay <- as.nifti(overlay, ...)
        }
    } else {
        if (is.character(overlays))
            overlays <- read.nifti(overlays)
        outs[[1]]$overlays.opts <- overlays.opts
        outs[[1]]$overlay <- as.nifti(overlays, ...)
    }
    
    # Threshold
    if (auto.thresh == TRUE) {
        if (is.null(thresh))
            thresh <- "overlay<=0"
        else
            warning("thresh is not NULL, so not using auto.thresh")
    }
    
    if (!is.null(thresh)) {
        for (o in 1:length(outs)) {
            overlay <- outs[[o]]$overlay
            thresh <- sprintf("tmp = %s", thresh)
            thresh <- eval(parse(text=thresh))
            outs[[o]]$overlay[thresh] <- NA
        }
    }
    
    # Binarize
    if (to.binarize) {
        for (o in 1:length(outs))
            outs[[o]]$overlay <- (outs[[o]]$overlay>0)*1
    }
    
    # Set the color range
    if (auto.colrange) {
        if (is.null(colrange))
            colrange <- range(sapply(outs, function(out) range(out$overlay, na.rm=T)))
        else
            warning("colrange is not NULL, so not using auto.colrange")
    } else if (is.null(colrange)) {
        stop("colrange cannot be NULL when auto.colrange is FALSE")
    }
    
    for (o in 1:length(outs)) {
        # 1. Set colrange
        outs[[o]]$overlays.opts$colrange <- colrange
        
        # 2. Get values that might not be shown
        overlay <- outs[[o]]$overlay
        values_notna <- !is.na(overlay)
        values_below <- overlay < colrange[1] & values_notna
        values_above <- overlay > colrange[2] & values_notna
        anyvalues_below <- any(values_below)
        anyvalues_above <- any(values_above)
        
        # 3. Check if any of these values exists
        if (anyvalues_below | anyvalues_above) {
            ## change them
            if (show.all.values) {
                warning("Values exist outside colrange...adjusting them")
                if (anyvalues_below)
                    outs[[o]]$overlay[values_below] <- colrange[1]
                if (anyvalues_above)
                    outs[[o]]$overlay[values_above] <- colrange[2]
            } else {
                warning("Values exist outside colrange...but NOT adjusting them")
            }
        }
    }
    
    # Return (basically list with each element having: col, colrange, and overlay)
    return(outs)
}

# Will return a specific slice in a 3D image
#' @nord
slice.nifti = function(image, orientation, slice) {
    # check that image is 3D
    if (length(dim(image)) != 3)
        stop("slice.nifti requires a 3D nifti image")
    
    switch(orientation,
        sagittal = image[slice,,],
        coronal = image[,slice,],
        axial = image[,,slice],
        stop("error orientation '", orientation, "' not recognized")
    )
}

#' @nord
prepare_slice_numbers <- function(slices, orientations, idim) {
    # 3 possibilities for slices
    ## 1. character => num of slices
    ## 2+3 list composed of either (can be intermixed, no problem)
    ### 2. negative number = absolute slice number
    ### 3. positive number that is 0-1 = relative slice number

    # TODO:
    # 1. have option for using center of mass for overlay as the slices
    # 2. have option to be able to use standard space units for slices (in this case would need to pass the header file and not just idim)
    
    get.idim <- function(idim, orientation) {
        switch(orientation,
            sagittal = idim[1],
            coronal = idim[2],
            axial = idim[3],
            stop("error orientation '", orientation, "' not recognized")
        )
    }
    
    # character (number of slices desired)
    if (is.character(slices)) {
        nslices <- as.integer(slices)
        if (length(orientations)!=1)
            stop("number of orientations given is larger than 1 (nslices)")
        idim <- get.idim(idim, orientations)
        slices <- seq(1, nslices, by=round(idim/nslices))*-1
    } 
    
    # orientation
    n <- length(slices)
    if (length(orientations) < n)
        orientations <- rep(orientations, length.out=n)
    else if (length(orientations) > n)
        stop("more orientations than necessary were given, number of orientation options should be length of number of slices or less")
    
    # convert any slices that are fractions to absolute #s
    # also check for bad slice #
    for (i in 1:n) {
        if (slices[i]>0 && slices[i]<1)
            slices[i] <- slices[i]*get.idim(idim, orientations[i])
        else if (slices[i]>get.idim(idim, orientations[i]))
            stop("slice #", slices[i], "is out of bounds in ", orientations[i], " orientation")
    }
    
    return(list(slices=slices, orientations=orientations))
}

# converts a slice between native and standard space
#' @nord
convert_slice <- function(header, ctype, slice, orientation) {
    tmat <- switch(ctype,
        native = header@qto.ijk,
        standard = header@qto.xz,
        stop("invalid ctype argument, must be native or standard")
    )
    
    new_slice <- switch(orientation,
        sagittal = c(slice,1) * tmat[1,c(1,4)],
        coronal = c(slice,1) * tmat[2,c(2,4)],
        axial = c(slice,1) * tmat[3,c(3,4)]
    )
    
    return(new_slice)
}

# plot the images
# TODO: add legend options and functions!!!
#' @nord
show.nimage <- function(p) {
    if (is.null(p$underlay) && is.null(p$underlay.opts))
        stop("show.nimage requires an underlay and underlay.opts to be specified")
    
    # 1. Defaults
    ## slice
    if (is.null(p$slice))
        p <- p + ni_slice()
    ## layout
    if (is.null(p$layout)) 
        p <- p + ni_layout()
    
    # 2. Layout
    nslices <- length(p$slice$slices)
    if (is.null(p$layout$nc))
        nc <- floor(sqrt(nslices))
    nr <- ceiling(nslices/nc)
    tmp <- c(1:nslices, rep(0, nc*nr - nslices))
    nf <- layout(matrix(tmp, nc, nr, byrow=TRUE))
    
    # 3. Check that overlays are all the same dimensions
    noverlays <- length(p$overlays)
    for (o in 2:noverlays) {
        if (!all( dim(p$overlays[[1]]) == dim(p$overlays[[2]]) ))
            stop("overlay dimensionsions do not match each other")
    }
    overlays <- p$overlays
    template_overlay <- overlays[[1]]
    
    # 4. Prepare slice numbers and orientations
    ## slice numbers related to underlay if no overlay
    if (noverlays == 0)
        idim <- dim(p$underlay)
    else
        idim <- dim(template_overlay)
    tmp <- prepare_slice_numbers(p$slice$slices, p$slice$orientations, idim)
    slices <- tmp$slices
    orientations <- tmp$orientations    
    
    # 5. Prepare margins and stuff
    if (nslices == 1) {
        par(mar=c(0,0.5,0,0.5))
        par(oma=c(0,0,0,0))
    } else {
        par(mar=c(1,1,1,1))
    }
    
    # 6. Loop through slices
    for (s in 1:nslices) {
        slice <- slices[s]
        orientation <- orientations[s]
        
        ## create underlay image
        ## if no overlays, then slice refers to underaly otherwise to overlays so need to convert
        if (noverlays == 0) {
            underlay_image <- slice.nifti(p$underlay, orientation, slice)
        } else {
            ## need to convert slice from overlay to underlay
            tmp_slice <- convert_slice(header(template_overlay), "standard", slice, orientation)
            tmp_slice <- convert_slice(header(p$underlay), "native", tmp_slice, orientation)
            underlay_image <- slice.nifti(p$underlay, orientation, round(tmp_slice))
        }
        
        ## plot underlay
        image(underlay_image, col=p$underlay.opts$col, axes=F, asp=ncol(underlay_image)/nrow(underlay_image))
        
        if (noverlays == 0)
            next
        
        ## plot overlays
        for (o in 1:noverlays) {
            overlay_image <- slice.nifti(overlays[[o]], orientation, slice)
            image(overlay_image, col=p$overlays.opts$col, add=T, axes=F, zlim=p$overlays.opts$colrange)
        }
        
        ## add slice numbers (in standard space)
        if (p$slice$show.numbers) {
            ## need to convert into standard space
            tmp_slice <- convert_slice(header(template_overlay), "standard", slice, orientation)
            txt = switch(orientation,
                sagittal = "X",
                axial = "Y",
                coronal = "Z"
            )
            par(oma=c(2,1,1,1))
            mtext(sprintf("%s = %i", txt, tmp_slice), 1, adj=1, font=2, cex=2) 
        }
    }    
}
