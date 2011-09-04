#include "connectir/connectir.h"

// Cs <- (t(As) %*% Bs)/(nrow(Bs)-1)
SEXP big_cor(SEXP As, SEXP Bs, SEXP Cs) {    
    try {        
        BM_TO_ARMA_INIT()
        BM_TO_ARMA_MULTIPLE(As, A)
        BM_TO_ARMA_MULTIPLE(Bs, B)
        double df = pMat->nrow() - 1;
        BM_TO_ARMA_MULTIPLE(Cs, C)
        
        C = (arma::trans(A) * B)/df;

        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

// Cs <- (As %*% t(Bs))/(nrow(As)-1)
SEXP big_tcor(SEXP As, SEXP Bs, SEXP Cs) {  
    try {        
        BM_TO_ARMA_INIT()
        BM_TO_ARMA_MULTIPLE(As, A)
        double df = pMat->nrow() - 1;
        BM_TO_ARMA_MULTIPLE(Bs, B)
        BM_TO_ARMA_MULTIPLE(Cs, C)
        
        C = (A * arma::trans(B))/df;

        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

// ztransform
SEXP big_ztransform(SEXP Cs) {
    try {
        BM_TO_ARMA_ONCE(Cs, C)
        atanh(C);
        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}
