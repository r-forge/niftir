#include "connectir/connectir.h"

// Cs <- (t(As) %*% Bs)/(nrow(Bs)-1)
SEXP big_cor(SEXP As, SEXP Bs, SEXP Cs, SEXP Cicol, SEXP Cncol) {    
    try {
        index_type icol = static_cast<index_type>(DOUBLE_DATA(Cicol)[0]);
        index_type ncol = static_cast<index_type>(DOUBLE_DATA(Cncol)[0]);
        
        BM_TO_ARMA_INIT()
        BM_TO_ARMA_MULTIPLE(As, A)
        BM_TO_ARMA_MULTIPLE(Bs, B)
        double df = pMat->nrow() - 1;
        
        bm = Cs;
        addr = bm.slot("address");
        pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
        if (pMat->matrix_type() != 8)
            ::Rf_error("Big Matrix must be of type double");
        offset = pMat->nrow() * pMat->col_offset();
        ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset + icol;
        arma::mat C(ptr_double, pMat->nrow(), ncol, false);
        
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
SEXP big_tcor(SEXP As, SEXP Bs, SEXP Cs, SEXP Cicol, SEXP Cncol) {  
    try {
        index_type icol = static_cast<index_type>(DOUBLE_DATA(Cicol)[0]);
        index_type ncol = static_cast<index_type>(DOUBLE_DATA(Cncol)[0]);
        
        BM_TO_ARMA_INIT()
        BM_TO_ARMA_MULTIPLE(As, A)
        BM_TO_ARMA_MULTIPLE(Bs, B)
        double df = pMat->nrow() - 1;
        
        bm = Cs;
        addr = bm.slot("address");
        pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
        if (pMat->matrix_type() != 8)
            ::Rf_error("Big Matrix must be of type double");
        offset = pMat->nrow() * pMat->col_offset();
        ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset + icol;
        arma::mat C(ptr_double, pMat->nrow(), ncol, false);
        
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
