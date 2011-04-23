#' @nord
.onLoad <- function(libname, pkgname) {
    library.dynam("bigextensions", pkgname, libname);
}

#.noGenerics <- TRUE           # This was a problem, not used.

#' @nord
.onUnload <- function(libpath) {
    library.dynam.unload("bigextensions", libpath);
}