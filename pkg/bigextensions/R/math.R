setMethod('tanh', 
    signature(x="big.matrix"),
    function(x) {
        .Call("BigTanhMain", x@address)
    }
)

setMethod('atanh', 
    signature(x="big.matrix"),
    function(x) {
        limit <- atanh(1-10^(-.Machine$double.exponent))
        .Call("BigAtanhMain", x@address, as.double(limit))
    }
)

setMethod('log', 
    signature(x="big.matrix"),
    function(x) {
        .Call("BigLogMain", x@address)
    }
)

setMethod('log10', 
    signature(x="big.matrix"),
    function(x) {
        .Call("BigLog10Main", x@address)
    }
)

setMethod('abs', 
    signature(x="big.matrix"),
    function(x) {
        .Call("BigAbsMain", x@address)
    }
)

setMethod('sqrt', 
    signature(x="big.matrix"),
    function(x) {
        .Call("BigSqrtMain", x@address)
    }
)


setGeneric("pow",
    function(x, exponent)
        standardGeneric('pow')
)

setMethod('pow', 
    signature(x="big.matrix"),
    function(x, exponent) {
        .Call("BigPowMain", x@address, as.double(exponent))
    }
)
