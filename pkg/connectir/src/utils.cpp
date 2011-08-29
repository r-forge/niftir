#include "connectir/connectir.h"

// Return row sum of big matrix (must be of type double)
SEXP bm_rowsum(SEXP Rbigmat) {
    SEXP res = R_NilValue;
    
    try {
        BM_TO_ARMA(Rbigmat, amat)
        using namespace arma;
        res = Rcpp::wrap(sum(amat, 1));
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return res;
}

// Return row mean of big matrix (must be of type double)
SEXP bm_rowmean(SEXP Rbigmat) {
    SEXP res = R_NilValue;
    
    try {
        BM_TO_ARMA(Rbigmat, amat)
        using namespace arma;
        res = Rcpp::wrap(mean(amat, 1));
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return res;
}
