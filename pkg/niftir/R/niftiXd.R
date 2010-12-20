#' @nord
setGeneric('is.niftiXd', 
    function(x) standardGeneric('is.niftiXd')
)

#' @nord
setMethod('is.niftiXd',
    signature(x='niftiXd'),
    function(x) return(TRUE)
)

#' @nord
setMethod('is.niftiXd',
    definition=function(x) return(FALSE)
)
