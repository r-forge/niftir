#include "connectir/connectir.h"

BigMatrix* sbm_to_bm(SEXP Sbm) {
    try {  
        Rcpp::RObject Rbm(Sbm);
        SEXP addr = Rbm.slot("address");
        BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
        return pMat;        
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return NULL;
}

void sbm_to_arma_xd(SEXP SbM, arma::mat& M) {
    try {
        BigMatrix *pMat = sbm_to_bm(SbM);
        if (pMat->matrix_type() != 8)
            Rf_error("big matrix must be of type double");
        if (pMat->row_offset() > 0 && pMat->ncol() > 1)
            Rf_error("big matrix cannot have a row offset and more than one column");
        
        index_type offset = pMat->nrow() * pMat->col_offset();
        double *ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
        
        arma::access::rw(M.mem) = ptr_double;
        arma::access::rw(M.n_rows) = pMat->nrow();
        arma::access::rw(M.n_cols) = pMat->ncol();
        arma::access::rw(M.n_elem) = pMat->nrow() * pMat->ncol();
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        Rf_error("c++ exception (unknown reason)");
    }
}

SEXP test_func(SEXP SbM) {
    try {
        arma::mat M;
        sbm_to_arma_xd(SbM, M);
        M = M * 2;
        return Rcpp::wrap( M );
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

void sub_sbm_to_arma_xd(SEXP SbM, arma::mat& M, SEXP SfirstCol, SEXP SlastCol) 
{
    try {
        index_type firstCol = static_cast<index_type>(DOUBLE_DATA(SfirstCol)[0] - 1);
        index_type lastCol = static_cast<index_type>(DOUBLE_DATA(SlastCol)[0] - 1);
        if (firstCol > lastCol)
            Rf_error("first column cannot be greater than the last column");
        index_type ncol = lastCol - firstCol + 1;
        
        BigMatrix *pMat = sbm_to_bm(SbM);
        if (pMat->matrix_type() != 8)
            Rf_error("big matrix must be of type double");
        if (pMat->row_offset() > 0 && pMat->ncol() > 1)
            Rf_error("big matrix cannot have a row offset and more than one column");
        if (lastCol > pMat->ncol())
            Rf_error("last column cannot be greater than # of columns in big matrix");
        
        index_type offset = pMat->nrow() * (pMat->col_offset() + firstCol);
        double *ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
        
        arma::access::rw(M.mem) = ptr_double;
        arma::access::rw(M.n_rows) = pMat->nrow();
        arma::access::rw(M.n_cols) = ncol;
        arma::access::rw(M.n_elem) = pMat->nrow() * ncol;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        Rf_error("c++ exception (unknown reason)");
    }
}

//
//SEXP test_func(SEXP Sbm) {
//    try {
//        BigMatrix *pMat = rbm_to_bm_xd(Sbm);
//        index_type offset = pMat->nrow() * pMat->col_offset();
//        double *ptr_double = reinterpret_cast<double*>(pMat->matrix()) + offset;
//        arma::mat amat(ptr_double, pMat->nrow(), pMat->ncol() - offset, false);
//        amat = amat * 2;
//        return R_NilValue;
//    } catch(std::exception &ex) {
//        forward_exception_to_r(ex);
//    } catch(...) {
//        ::Rf_error("c++ exception (unknown reason)");
//    }
//    return R_NilValue;
//}
