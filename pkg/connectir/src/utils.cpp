#include "connectir/connectir.h"

//BigMatrix* rbm_to_bm_xd(SEXP Sbm) {
//    try {  
//        Rcpp::RObject Rbm(Sbm);
//        SEXP addr = Rbm.slot("address");
//        BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
//    
//        if (pMat->matrix_type() != 8)
//            Rf_error("Big Matrix must be of type double");
//        
//        return pMat;        
//    } catch(std::exception &ex) {
//        forward_exception_to_r(ex);
//    } catch(...) {
//        ::Rf_error("c++ exception (unknown reason)");
//    }
//}
//
//SEXP test_func(SEXP Sbm) {
//    try {
//        BigMatrix *pMat = rbm_to_bm_xd(Sbm);
//        index_type offset = pMat->nrow() * pMat->col_offset();
//        double *ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
//        arma::mat amat(ptr_double, pMat->nrow(), pMat->ncol() - offset, false);
//        amat * 2;
//        return R_NilValue;
//    } catch(std::exception &ex) {
//        forward_exception_to_r(ex);
//    } catch(...) {
//        ::Rf_error("c++ exception (unknown reason)");
//    }
//    return R_NilValue;
//}

// Return row sum of big matrix (must be of type double)
// Y = a*X + b
SEXP big_add_multiply_scalar(SEXP SX, SEXP SY,  
                             SEXP Sa, SEXP Sb, 
                             SEXP SX_firstCol, SEXP SX_lastCol, 
                             SEXP SY_firstCol, SEXP SY_lastCol) 
{
    try {        
        double a = DOUBLE_DATA(Sa)[0];
        double b = DOUBLE_DATA(Sb)[0];
        
        BM_COL_INIT(SX_firstCol, X_firstCol)
        BM_COL_INIT(SX_lastCol, X_lastCol)
        BM_COL_INIT(SY_firstCol, Y_firstCol)
        BM_COL_INIT(SY_lastCol, Y_lastCol)
        
        BM_TO_ARMA_INIT()
        SUB_BM_TO_ARMA_MULTIPLE(SX, X, X_firstCol, X_lastCol)
        SUB_BM_TO_ARMA_MULTIPLE(SY, Y, Y_firstCol, Y_lastCol)
        
        Y = a*X + b;
        
        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

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
