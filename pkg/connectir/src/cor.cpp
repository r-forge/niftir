#include "connectir/connectir.h"

// In R:
// library(connectir)
// x <- as.big.matrix(matrix(1:20, 5, 4), type="double")
// y <- .Call("test_sub_matrix", x, as.double(2), as.double(4))
// all.equal(x[,2:4], y)
SEXP test_sub_matrix(SEXP As, SEXP As_firstCol, SEXP As_lastCol) {
    try {
        BM_COL_INIT(As_firstCol, A_firstCol)
        BM_COL_INIT(As_lastCol, A_lastCol)
    
        BM_TO_ARMA_INIT()
        
        SUB_BM_TO_ARMA_MULTIPLE(As, A, A_firstCol, A_lastCol)
        
        return Rcpp::wrap( A );
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

// Cs <- (t(As) %*% Bs)/(nrow(Bs)-1)
SEXP big_cor(SEXP As, SEXP Bs, SEXP Cs, 
             SEXP As_firstCol, SEXP As_lastCol, 
             SEXP Bs_firstCol, SEXP Bs_lastCol, 
             SEXP Cs_firstCol, SEXP Cs_lastCol) 
{    
    try {
        BM_COL_INIT(As_firstCol, A_firstCol)
        BM_COL_INIT(As_lastCol, A_lastCol)
        BM_COL_INIT(Bs_firstCol, B_firstCol)
        BM_COL_INIT(Bs_lastCol, B_lastCol)
        BM_COL_INIT(Cs_firstCol, C_firstCol)
        BM_COL_INIT(Cs_lastCol, C_lastCol)
        
        BM_TO_ARMA_INIT()
        SUB_BM_TO_ARMA_MULTIPLE(As, A, A_firstCol, A_lastCol)
        SUB_BM_TO_ARMA_MULTIPLE(Bs, B, B_firstCol, B_lastCol)
        double df = pMat->nrow() - 1;
        SUB_BM_TO_ARMA_MULTIPLE(Cs, C, C_firstCol, C_lastCol)
        
        // TODO: make this an explicit user option!
        if (ncol == 1 && B.n_cols != 1 && C.n_rows == (A.n_cols*B.n_cols))
            C.reshape(A.n_cols, B.n_cols);
        
        if (A.n_cols != C.n_rows)
            Rf_error("incorrect number of rows in output matrix");
        if (B.n_cols != C.n_cols)
            Rf_error("incorrect number of columns in output matrix");
        
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
SEXP big_tcor(SEXP As, SEXP Bs, SEXP Cs, 
             SEXP As_firstCol, SEXP As_lastCol, 
             SEXP Bs_firstCol, SEXP Bs_lastCol, 
             SEXP Cs_firstCol, SEXP Cs_lastCol) 
{    
    try {
        BM_COL_INIT(As_firstCol, A_firstCol)
        BM_COL_INIT(As_lastCol, A_lastCol)
        BM_COL_INIT(Bs_firstCol, B_firstCol)
        BM_COL_INIT(Bs_lastCol, B_lastCol)
        BM_COL_INIT(Cs_firstCol, C_firstCol)
        BM_COL_INIT(Cs_lastCol, C_lastCol)
        
        BM_TO_ARMA_INIT()
        SUB_BM_TO_ARMA_MULTIPLE(As, A, A_firstCol, A_lastCol)
        SUB_BM_TO_ARMA_MULTIPLE(Bs, B, B_firstCol, B_lastCol)
        double df = ncol - 1;
        SUB_BM_TO_ARMA_MULTIPLE(Cs, C, C_firstCol, C_lastCol)
        
        // TODO: make this an explicit user option!
        if (ncol == 1 && B.n_rows != 1 && C.n_rows == (A.n_rows*B.n_rows))
            C.reshape(A.n_rows, B.n_rows);
        
        if (A.n_rows != C.n_rows)
            Rf_error("incorrect number of rows in output matrix");
        if (B.n_rows != C.n_cols)
            Rf_error("incorrect number of columns in output matrix");
        
        C = (A * arma::trans(B))/df;
        
        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

