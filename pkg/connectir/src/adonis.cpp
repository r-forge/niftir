#include "connectir/connectir.h"

SEXP mdmr_nmat_to_pmat(SEXP SNmat, SEXP Snperms)
{
    try {
        arma::mat Nmat(1,1);
        const double* old_nptr = sbm_to_arma_xd(SNmat, Nmat);
        double nperms = DOUBLE_DATA(Snperms)[0];
        Nmat /= nperms;
        free_arma(Nmat, old_nptr);
        return SNmat;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

SEXP mdmr_worker(SEXP SfirstVox, SEXP SlastVox, SEXP SGmat, SEXP SFperms, 
                 SEXP SNmat, 
                 SEXP SH2mats, SEXP SIHmat, 
                 SEXP SdfRes, SEXP SdfExp)
{
    try {
        Rcpp::List tmp(SH2mats);
        index_type nterms = static_cast<index_type>(tmp.size());
        
        arma::mat Gmat(1,1);
        const double* old_gptr = sub_sbm_to_arma_xd(SGmat, Gmat, SfirstVox, SlastVox);
        index_type nvoxs = static_cast<index_type>(Gmat.n_cols);
        
        SEXP SFmat;
        arma::mat Fmat(1,1); const double* old_fptr = Fmat.memptr();
        
        arma::mat Nmat(1,1);
        const double* old_nptr = sub_sbm_to_arma_xd(SNmat, Nmat, SfirstVox, SlastVox);
        arma::rowvec realFs; double* Nvec;
        
        SEXP SH2mat;
        arma::mat H2mat(1,1); const double* old_h2ptr = H2mat.memptr();
        
        arma::mat IHmat(1,1);
        const double* old_ihptr = sbm_to_arma_xd(SIHmat, IHmat);
        
        double dfRes = DOUBLE_DATA(SdfRes)[0];
        Rcpp::NumericVector RdfExp(SdfExp);
        arma::vec dfExp(RdfExp.begin(), RdfExp.size(), false);
        
        arma::mat ExplainedVariance;
        arma::mat ErrorVariance = arma::trans(IHmat) * Gmat;
        for (index_type i=0; i < nterms; ++i)
        {
            PROTECT(SH2mat = VECTOR_ELT(SH2mats, i));
            sbm_to_arma_xd(SH2mat, H2mat);
            UNPROTECT(1);

            PROTECT(SFmat = VECTOR_ELT(SFperms, i));
            sub_sbm_to_arma_xd(SFmat, Fmat, SfirstVox, SlastVox);
            UNPROTECT(1);
            
            ExplainedVariance = arma::trans(H2mat) * Gmat;
            
            // ExplainedVariance/ErrorVariance
            Fmat = (ExplainedVariance/ErrorVariance) * (dfRes/dfExp(i));
            
            // # of observations greater than original
            realFs = Fmat.row(0);
            Nvec = const_cast<double *>(Nmat.colptr(i));
            for (i = 0; i < nvoxs; ++i)
            {
                Nvec[i] += arma::as_scalar(arma::sum(realFs(i) >= Fmat.col(i)));
            }
        }
        
        free_arma(Gmat, old_gptr);
        free_arma(Nmat, old_nptr);
        free_arma(H2mat, old_h2ptr);
        free_arma(IHmat, old_ihptr);
        free_arma(Fmat, old_fptr);
        
        return SFperms;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

template<typename CType, typename BMAccessorType>
SEXP ComputePvals(BigMatrix *inMat, BigMatrix *outMat, double colnum) {
    BMAccessorType om( *outMat );
    CType *outCol;
    outCol = om[static_cast<index_type>(colnum)-1];
    
    BMAccessorType im( *inMat );    
    index_type ncols = inMat->ncol();
    index_type nrows = inMat->nrow();
    index_type i=0;
    index_type j=0;
    CType *inCol;
    
    for (i=0; i < ncols; ++i) {
        inCol = im[i];
        outCol[i] = 1;
        for (j=1; j < nrows; ++j) {
            if (inCol[j] > inCol[0])
                outCol[i] += 1;
        }
        outCol[i] = outCol[i]/nrows;
    }
    
    return R_NilValue;
}


extern "C" {

SEXP ComputePvalsMain(SEXP Rinmat, SEXP Routmat, SEXP Routcol) {
    BigMatrix *inMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(Rinmat));
    BigMatrix *outMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(Routmat));
    double outCol = NUMERIC_DATA(Routcol)[0];
    
    if (inMat->separated_columns() != outMat->separated_columns())
        Rf_error("all big matrices are not the same column separated type");
    if (inMat->matrix_type() != outMat->matrix_type())
        Rf_error("all big matrices are not the same matrix type");
    if (inMat->ncol() != outMat->nrow())
        Rf_error("inMat # of cols must be the same as outMat # of rows");
    
    CALL_BIGFUNCTION_ARGS_THREE(ComputePvals, inMat, outMat, outCol)
    return(ret);
}

} // end extern "C"
