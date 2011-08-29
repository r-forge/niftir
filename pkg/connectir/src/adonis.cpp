#include "connectir/connectir.h"

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
