#include "connectir/connectir.h"

//arma::mat bm_to_arma_d(SEXP Sbm) {
//    Rcpp::RObject Rbm(Sbm);
//    SEXP addr = Rbm.slot("address");
//    BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
//    
//    if (pMat->matrix_type() != 8)
//        Rf_error("Big Matrix must be of type double");
//    
//    index_type offset = pMat->nrow() * pMat->col_offset();
//    double *ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
//    
//    arma::mat amat(ptr_double, pMat->nrow(), pMat->ncol(), false);
//    return amat;
//}

// Return row sum of big matrix (must be of type double)
SEXP bm_rowsum(SEXP Rbigmat) {
    SEXP res = R_NilValue;
    
    try {
        BM_TO_ARMA_ONCE(Rbigmat, amat)
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
        BM_TO_ARMA_ONCE(Rbigmat, amat)
        using namespace arma;
        res = Rcpp::wrap(mean(amat, 1));
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return res;
}
