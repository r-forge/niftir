vbca <- function(bigmat, cols, ztransform=FALSE, outmat=NULL, ...) {
    A <- deepcopy(bigmat, cols=cols)
#    A <- sub.big.matrix(bigmat, firstCol=firstCol, lastCol=lastCol)
	B <- bigmat
	if (is.null(outmat))
	    C <- big.matrix(ncol(A), ncol(B), init=0, type="double", ...)
	else if (ncol(A)==nrow(outmat) && ncol(B)==ncol(outmat))
	    C <- outmat
	else
	    stop("dimensions of outmat are out of whack")
	ALPHA <- 1/(nrow(B) - 1)
	dgemm(C=C, A=A, B=B, TRANSA='t', ALPHA=ALPHA)
	if (ztransform)
	    atanh(C)
	invisible(C)
}

vbca_batch <- function(subs.bigmats, cols, ztransform=FALSE, ...) {
    nsubjects <- length(subs.bigmats)
    lapply(1:nsubjects, function(i) vbca(subs.bigmats[[i]], cols, ztransform, ...))
}



## NEW CODE USING ARMADILLO

library(RcppArmadillo)
library(inline)

plugin_bigmemory <- function() {
  l <- getPlugin("RcppArmadillo")
  
  l$includes <- paste(l$includes, '
#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/bigmemoryDefines.h"
#include "bigmemory/isna.hpp"
')
  
  l$LinkingTo <- c("bigmemory", l$LinkingTo)
  
  l$Depends <- c("bigmemory", l$Depends)  
  
  return(l)
}
registerPlugin("bigmemory", plugin_bigmemory)

# c <- (t(A) %*% B)/(nrow(B)-1)
cpp_cor <- cxxfunction( signature(As="object", Bs="object", Cs="object"), 
'
  try{
    using namespace Rcpp;
    
    SEXP addr; BigMatrix *pMat; index_type offset; double *ptr_double;
    
    RObject Abm(As);
    addr = Abm.slot("address");
    pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
    offset = pMat->nrow() * pMat->col_offset();
    ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
    arma::mat A(ptr_double, pMat->nrow(), pMat->ncol(), false);
    
    RObject Bbm(Bs);
    addr = Bbm.slot("address");
    pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
    offset = pMat->nrow() * pMat->col_offset();
    ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
    arma::mat B(ptr_double, pMat->nrow(), pMat->ncol(), false);
    double df = pMat->nrow() - 1;
    
    RObject Cbm(Cs);
    addr = Cbm.slot("address");
    pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
    offset = pMat->nrow() * pMat->col_offset();
    ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
    arma::mat C(ptr_double, pMat->nrow(), pMat->ncol(), false);
    
    C = (arma::trans(A) * B)/df;
    
    return R_NilValue;
  } catch( std::exception &ex ) {
  forward_exception_to_r( ex );
  } catch(...) { 
  ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
', plugin = "bigmemory")

# c <- (A %*% t(B))/(nrow(A)-1)
cpp_tcor <- cxxfunction( signature(As="object", Bs="object", Cs="object"), 
'
  try{
    using namespace Rcpp;
    
    SEXP addr; BigMatrix *pMat; index_type offset; double *ptr_double;
    
    RObject Abm(As);
    addr = Abm.slot("address");
    pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
    offset = pMat->nrow() * pMat->col_offset();
    ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
    arma::mat A(ptr_double, pMat->nrow(), pMat->ncol(), false);
    
    RObject Bbm(Bs);
    addr = Bbm.slot("address");
    pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
    offset = pMat->nrow() * pMat->col_offset();
    ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
    arma::mat B(ptr_double, pMat->nrow(), pMat->ncol(), false);
    double df = pMat->nrow() - 1;
    
    RObject Cbm(Cs);
    addr = Cbm.slot("address");
    pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
    offset = pMat->nrow() * pMat->col_offset();
    ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
    arma::mat C(ptr_double, pMat->nrow(), pMat->ncol(), false);
    
    C = (A * arma::trans(B))/df;
    
    return R_NilValue;
  } catch( std::exception &ex ) {
  forward_exception_to_r( ex );
  } catch(...) { 
  ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
', plugin = "bigmemory")

cpp_ztransform <- cxxfunction( signature(Cs="object"), 
'
  try{
    using namespace Rcpp;
    
    SEXP addr; BigMatrix *pMat; index_type offset; double *ptr_double;
    
    RObject Cbm(Cs);
    addr = Cbm.slot("address");
    pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
    offset = pMat->nrow() * pMat->col_offset();
    ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
    arma::mat C(ptr_double, pMat->nrow(), pMat->ncol(), false);
    
    atanh(C);
    
    return R_NilValue;
  } catch( std::exception &ex ) {
  forward_exception_to_r( ex );
  } catch(...) { 
  ::Rf_error( "c++ exception (unknown reason)" ); 
  }
  return R_NilValue; // -Wall
', plugin = "bigmemory")


vbca2 <- function(bigmat, cols, ztransform=FALSE, outmat=NULL, ...) {
        A <- deepcopy(bigmat, cols=cols)
    #    A <- sub.big.matrix(bigmat, firstCol=firstCol, lastCol=lastCol)
    	B <- bigmat
    	if (is.null(outmat))
    	    C <- big.matrix(ncol(A), ncol(B), init=0, type="double", ...)
    	else if (ncol(A)==nrow(outmat) && ncol(B)==ncol(outmat))
    	    C <- outmat
    	else
    	    stop("dimensions of outmat are out of whack")
    	
    	cpp_cor(A, B, C)
    	if (ztransform)
    	    atanh(C)
    	invisible(C)
}

vbca_batch2 <- function(subs.bigmats, cols, ztransform=FALSE, ...) {
    nsubjects <- length(subs.bigmats)
    lapply(subs.bigmats, function(x) vbca2(x, cols, ztransform, ...))
}
