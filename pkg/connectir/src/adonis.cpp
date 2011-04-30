#include "connectir/connectir.h"

template<typename CType, typename BMAccessorType>
void ComputePvals(BigMatrix *pMat, double *pRet) {
    BMAccessorType m( *pMat );
    
    index_type ncols = pMat->ncol();
    index_type nrows = pMat->nrow();
    index_type i=0;
    index_type j=0;
    CType *pCol;
    
    for (i=0; i < ncols; ++i) {
        pCol = m[i];
        pRet[i] = 1;
        for (j=1; j < nrows; ++j) {
            if (pCol[j] > pCol[0])
                pRet[i] += 1;            
        }
        pRet[i] = pRet[i]/nrows;
    }
    
    return;
}


extern "C" {

SEXP ComputePvalsMain(SEXP bigaddr) {
    BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(bigaddr));
    index_type nCols = pMat->ncol();
    
    SEXP ret = R_NilValue;
    ret = PROTECT(NEW_NUMERIC(nCols));
    double *pRet = NUMERIC_DATA(ret);
    
    if (pMat->separated_columns()) { 
        switch (pMat->matrix_type()) {
            case 1:
                ComputePvals<char, SepMatrixAccessor<char> >(
                    pMat, pRet);
                break;
            case 2:
                ComputePvals<short, SepMatrixAccessor<short> >(
                    pMat, pRet);
                break;
            case 4:
                ComputePvals<int, SepMatrixAccessor<int> >(
                    pMat, pRet);
                break;
            case 8:
                ComputePvals<double, SepMatrixAccessor<double> >(
                    pMat, pRet);
                break;
        }
    }
    else {
        switch (pMat->matrix_type()) {
            case 1:
                ComputePvals<char, MatrixAccessor<char> >(
                    pMat, pRet);
                break;
            case 2:
                ComputePvals<short, MatrixAccessor<short> >(
                    pMat, pRet);
                break;
            case 4:
                ComputePvals<int, MatrixAccessor<int> >(
                    pMat, pRet);
                break;
            case 8:
                ComputePvals<double, MatrixAccessor<double> >(
                    pMat, pRet);
                break;
        }
    }
    
    UNPROTECT(1);
    return(ret);
}

} // end extern "C"
