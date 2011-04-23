test.scale <- function() {
    xorig <- matrix(rnorm(50*4), 50, 4)
    x1 <- scale(xorig)
    x2 <- as.big.matrix(xorig)
    .Call("BigScale2Main", x2@address, as.double(1:nrow(x2)), as.double(1:ncol(x2)))
    return(list(xorig, x1, x2))
}

test.scale.time <- function() {
    xorig <- matrix(rnorm(5000*400), 5000, 400)
    system.time(x1 <- scale(xorig))
    x2 <- as.big.matrix(xorig)
    system.time(.Call("BigScaleMain", x2@address, as.double(1:nrow(x2)), as.double(1:ncol(x2))))
    x2 <- as.big.matrix(xorig)
    system.time(.Call("BigScale2Main", x2@address, as.double(1:nrow(x2)), as.double(1:ncol(x2))))
    x2 <- as.big.matrix(xorig)
    system.time(scale(x2, to.copy=FALSE))
}

test.transpose <- function() {
    xmat <- matrix(rnorm(200), 20, 10)
    t.xmat <- t(xmat)
    xbigmat <- as.big.matrix(xmat)
    t.xbigmat <- big.matrix(ncol(xbigmat), nrow(xbigmat))
    .Call("BigTransposeMain", xbigmat@address, t.xbigmat@address)
    list(all.equal(as.vector(t.xmat), as.vector(t.xbigmat[,])), t.xmat, t.xbigmat)
}
