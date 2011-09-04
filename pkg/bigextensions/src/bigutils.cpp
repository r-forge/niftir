#ifndef BIGEXTENSIONS_UTILS
#define BIGEXTENSIONS_UTILS

#include "bigextensions/bigex.h"

template<typename CType, typename BMAccessorType>
void BigTranspose(BigMatrix *pOrigMat, BigMatrix *pNewMat) {
    BMAccessorType origmat( *pOrigMat );
    BMAccessorType newmat( *pNewMat );
    
    index_type numRows = pOrigMat->nrow();
    index_type numCols = pOrigMat->ncol();
    
    if (numRows != pNewMat->ncol())
        error("length of row indices does not equal the # of cols in the new matrix");
    if (numCols != pNewMat->nrow())
        error("length of col indices does not equal the # of rows in the new matrix");
    
    index_type i=0;
    index_type j=0;
    for (i = 0; i < numCols; ++i)
        for (j = 0; j < numRows; ++j)
            newmat[j][i] = origmat[i][j];
    
    return;
}

template<typename CType, typename BMAccessorType>
void GetDiag(BigMatrix *pMat, double *pValues, index_type n) {
    BMAccessorType m( *pMat );
    for (index_type i = 0; i < n; ++i)
        pValues[i] = (double)m[i][i];
    return;
}

template<typename CType, typename BMAccessorType>
void SetDiag(BigMatrix *pMat, double *pValues, index_type n) {
    BMAccessorType m( *pMat );
    for (index_type i = 0; i < n; ++i)
        m[i][i] = (CType)pValues[i];
    return;
}

template<typename CType, typename BMAccessorType>
void BigDeepCopy(BigMatrix *pOrigMat, BigMatrix *pNewMat, SEXP rowIndices, SEXP colIndices) {
    BMAccessorType origmat( *pOrigMat );
    BMAccessorType newmat( *pNewMat );
    
    double *pRows = NUMERIC_DATA(rowIndices);
    double *pCols = NUMERIC_DATA(colIndices);
    index_type numRows = GET_LENGTH(rowIndices);
    index_type numCols = GET_LENGTH(colIndices);
    
    if (numRows != pNewMat->nrow())
        error("length of row indices does not equal the # of rows in the new matrix");
    if (numCols != pNewMat->ncol())
        error("length of column indices does not equal the # of columns in the new matrix");
    
    index_type i=0;
    index_type j=0;
    CType *pOrigColumn;
    CType *pNewColumn;
    
    for (i = 0; i < numCols; ++i) {
        pOrigColumn = origmat[static_cast<index_type>(pCols[i])-1];
        pNewColumn = newmat[i];
        for (j = 0; j < numRows; ++j)
            pNewColumn[j] = pOrigColumn[static_cast<index_type>(pRows[j])-1];
    }
    
    return;
}

template<typename CType, typename BMAccessorType>
void BigScale(BigMatrix *pMat, SEXP rowIndices, SEXP colIndices) {
    BMAccessorType m( *pMat );
    
    double *pRows = NUMERIC_DATA(rowIndices);
    double *pCols = NUMERIC_DATA(colIndices);
    index_type numRows = GET_LENGTH(rowIndices);
    index_type numCols = GET_LENGTH(colIndices);
    
    index_type i=0;
    index_type j=0;
    index_type jj=0;
    CType *pColumn;
    
    LDOUBLE x = 0;
    LDOUBLE delta = 0;
    LDOUBLE mean = 0;
    LDOUBLE M2 = 0;
    LDOUBLE stdev = 0;
    LDOUBLE scaled_x;
    
    for (i = 0; i < numCols; ++i) {
        pColumn = m[static_cast<index_type>(pCols[i])-1];
        
        // First pass to get mean and sd
        delta = 0;
        mean = 0;
        M2 = 0;
        stdev = 0;
        for (j = 0; j < numRows; ++j) {
            // todo: add checking for NaN...but shouldn't really have any!
            // maybe can also pass the exact list of voxs to loop through!
            // if (!ISNAN(pColumn[curj]))
            // NA_REAL
            x = static_cast<LDOUBLE>(pColumn[static_cast<index_type>(pRows[j])-1]);
            delta = x - mean;
            mean = mean + delta/static_cast<LDOUBLE>(j+1);
            M2 = M2 + delta*(x - mean);
        }
        stdev = sqrt(M2/(static_cast<LDOUBLE>(numRows-1)));
        
        //printf("mean: %f; stdev: %f\n", mean, stdev);
        
        // Second pass to scale
        for (j = 0; j < numRows; ++j) {
            jj = static_cast<index_type>(pRows[j]-1);
            scaled_x = (static_cast<LDOUBLE>(pColumn[jj])-mean)/stdev;
            pColumn[jj] = static_cast<CType>(scaled_x);
        }
    }
    
    return;
}

template<typename CType, typename BMAccessorType>
void BigTransScale(BigMatrix *pMat, SEXP rowIndices, SEXP colIndices) {
    BMAccessorType m( *pMat );
    
    double *pRows = NUMERIC_DATA(rowIndices);
    double *pCols = NUMERIC_DATA(colIndices);
    index_type numRows = GET_LENGTH(rowIndices);
    index_type numCols = GET_LENGTH(colIndices);
    
    index_type i=0;
    index_type j=0;
    index_type ci=0;
    index_type ri=0;
    CType *pColumn;
    
    LDOUBLE x = 0;
    LDOUBLE delta = 0;
    LDOUBLE mean = 0;
    LDOUBLE M2 = 0;
    LDOUBLE stdev = 0;
    LDOUBLE scaled_x;
    
    for (j = 0; j < numRows; ++j) {
        ri = static_cast<index_type>(pRows[j])-1;
        
        // First pass to get mean and sd
        delta = mean = M2 = stdev = 0;
        for (i = 0; i < numCols; ++i) {
            ci = static_cast<index_type>(pCols[i])-1;            
            x = static_cast<LDOUBLE>(m[ci][ri]);
            delta = x - mean;
            mean = mean + delta/static_cast<LDOUBLE>(i+1);
            M2 = M2 + delta*(x - mean);
        }
        stdev = sqrt(M2/(static_cast<LDOUBLE>(numCols-1)));
        
        // Second pass to scale
        for (i = 0; i < numCols; ++i) {
            ci = static_cast<index_type>(pCols[i]-1);
            scaled_x = (static_cast<LDOUBLE>(m[ci][ri])-mean)/stdev;
            m[ci][ri] = static_cast<CType>(scaled_x);
        }
        
    }
    
    return;
}

extern "C" {
    
    SEXP BigTransposeMain(SEXP origAddr, SEXP newAddr) {

        BigMatrix *pOrigMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(origAddr));
        BigMatrix *pNewMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(newAddr));

        if (pOrigMat->separated_columns() != pNewMat->separated_columns() | 
            pOrigMat->matrix_type() != pNewMat->matrix_type())
            error("original and new big matrix are not the same type");

        if (pOrigMat->separated_columns()) {
            switch (pOrigMat->matrix_type()) {
                case 1:
                    BigTranspose<char, SepMatrixAccessor<char> >( 
                    pOrigMat, pNewMat);
                    break;
                case 2:
                    BigTranspose<short, SepMatrixAccessor<short> >( 
                    pOrigMat, pNewMat);
                    break;
                case 4:
                    BigTranspose<int, SepMatrixAccessor<int> >( 
                    pOrigMat, pNewMat);
                    break;
                case 8:
                    BigTranspose<double, SepMatrixAccessor<double> >( 
                    pOrigMat, pNewMat);
                    break;
            }
        }
        else {
            switch (pOrigMat->matrix_type()) {
                case 1:
                    BigTranspose<char, MatrixAccessor<char> >( 
                    pOrigMat, pNewMat);
                    break;
                case 2:
                    BigTranspose<short, MatrixAccessor<short> >( 
                    pOrigMat, pNewMat);
                    break;
                case 4:
                    BigTranspose<int, MatrixAccessor<int> >( 
                    pOrigMat, pNewMat);
                    break;
                case 8:
                    BigTranspose<double, MatrixAccessor<double> >( 
                    pOrigMat, pNewMat);
                    break;
            }
        }

        return R_NilValue;
    }
    
    SEXP GetDiagMain(SEXP addr, SEXP Rn) {

        BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
        
        index_type n = NUMERIC_DATA(Rn)[0];
        
        SEXP ret = R_NilValue;
        ret = PROTECT(NEW_NUMERIC(n));
        double *pRet = NUMERIC_DATA(ret);

        if (pMat->separated_columns()) {
            switch (pMat->matrix_type()) {
                case 1:
                    GetDiag<char, SepMatrixAccessor<char> >( 
                    pMat, pRet, n);
                    break;
                case 2:
                    GetDiag<short, SepMatrixAccessor<short> >( 
                    pMat, pRet, n);
                    break;
                case 4:
                    GetDiag<int, SepMatrixAccessor<int> >( 
                    pMat, pRet, n);
                    break;
                case 8:
                    GetDiag<double, SepMatrixAccessor<double> >( 
                    pMat, pRet, n);
                    break;
            }
        }
        else {
            switch (pMat->matrix_type()) {
                case 1:
                    GetDiag<char, MatrixAccessor<char> >( 
                    pMat, pRet, n);
                    break;
                case 2:
                    GetDiag<short, MatrixAccessor<short> >( 
                    pMat, pRet, n);
                    break;
                case 4:
                    GetDiag<int, MatrixAccessor<int> >( 
                    pMat, pRet, n);
                    break;
                case 8:
                    GetDiag<double, MatrixAccessor<double> >( 
                    pMat, pRet, n);
                    break;
            }
        }

        UNPROTECT(1);
        return(ret);
    }
    
    SEXP SetDiagMain(SEXP addr, SEXP values) {

        BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
        
        double *pValues = NUMERIC_DATA(values);
        index_type n = GET_LENGTH(values);

        if (pMat->separated_columns()) {
            switch (pMat->matrix_type()) {
                case 1:
                    SetDiag<char, SepMatrixAccessor<char> >( 
                    pMat, pValues, n);
                    break;
                case 2:
                    SetDiag<short, SepMatrixAccessor<short> >( 
                    pMat, pValues, n);
                    break;
                case 4:
                    SetDiag<int, SepMatrixAccessor<int> >( 
                    pMat, pValues, n);
                    break;
                case 8:
                    SetDiag<double, SepMatrixAccessor<double> >( 
                    pMat, pValues, n);
                    break;
            }
        }
        else {
            switch (pMat->matrix_type()) {
                case 1:
                    SetDiag<char, MatrixAccessor<char> >( 
                    pMat, pValues, n);
                    break;
                case 2:
                    SetDiag<short, MatrixAccessor<short> >( 
                    pMat, pValues, n);
                    break;
                case 4:
                    SetDiag<int, MatrixAccessor<int> >( 
                    pMat, pValues, n);
                    break;
                case 8:
                    SetDiag<double, MatrixAccessor<double> >( 
                    pMat, pValues, n);
                    break;
            }
        }
        
        return R_NilValue;
    }
    
    SEXP BigDeepCopyMain(SEXP origAddr, SEXP newAddr, SEXP rowIndices, SEXP colIndices) {

        BigMatrix *pOrigMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(origAddr));
        BigMatrix *pNewMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(newAddr));

        if (pOrigMat->separated_columns() != pNewMat->separated_columns() | 
            pOrigMat->matrix_type() != pNewMat->matrix_type())
            error("original and new big matrix are not the same type");

        if (pOrigMat->separated_columns()) {
            switch (pOrigMat->matrix_type()) {
                case 1:
                    BigDeepCopy<char, SepMatrixAccessor<char> >( 
                    pOrigMat, pNewMat, rowIndices, colIndices);
                    break;
                case 2:
                    BigDeepCopy<short, SepMatrixAccessor<short> >( 
                    pOrigMat, pNewMat, rowIndices, colIndices);
                    break;
                case 4:
                    BigDeepCopy<int, SepMatrixAccessor<int> >( 
                    pOrigMat, pNewMat, rowIndices, colIndices);
                    break;
                case 8:
                    BigDeepCopy<double, SepMatrixAccessor<double> >( 
                    pOrigMat, pNewMat, rowIndices, colIndices);
                    break;
            }
        }
        else {
            switch (pOrigMat->matrix_type()) {
                case 1:
                    BigDeepCopy<char, MatrixAccessor<char> >( 
                    pOrigMat, pNewMat, rowIndices, colIndices);
                    break;
                case 2:
                    BigDeepCopy<short, MatrixAccessor<short> >( 
                    pOrigMat, pNewMat, rowIndices, colIndices);
                    break;
                case 4:
                    BigDeepCopy<int, MatrixAccessor<int> >( 
                    pOrigMat, pNewMat, rowIndices, colIndices);
                    break;
                case 8:
                    BigDeepCopy<double, MatrixAccessor<double> >( 
                    pOrigMat, pNewMat, rowIndices, colIndices);
                    break;
            }
        }

        return R_NilValue;
    }
    
    SEXP BigScaleMain(SEXP addr, SEXP rowIndices, SEXP colIndices) {

        BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
        if (pMat->separated_columns()) {
            switch (pMat->matrix_type()) {
                case 1:
                    BigScale<char, SepMatrixAccessor<char> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 2:
                    BigScale<short, SepMatrixAccessor<short> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 4:
                    BigScale<int, SepMatrixAccessor<int> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 8:
                    BigScale<double, SepMatrixAccessor<double> >( 
                    pMat, rowIndices, colIndices);
                    break;
            }
        }
        else {
            switch (pMat->matrix_type()) {
                case 1:
                    BigScale<char, MatrixAccessor<char> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 2:
                    BigScale<short, MatrixAccessor<short> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 4:
                    BigScale<int, MatrixAccessor<int> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 8:
                    BigScale<double, MatrixAccessor<double> >( 
                    pMat, rowIndices, colIndices);
                    break;
            }
        }
        
        return R_NilValue;
    }
    
    SEXP BigTransScaleMain(SEXP addr, SEXP rowIndices, SEXP colIndices) {
        BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));
        if (pMat->separated_columns()) {
            switch (pMat->matrix_type()) {
                case 1:
                    BigTransScale<char, SepMatrixAccessor<char> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 2:
                    BigTransScale<short, SepMatrixAccessor<short> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 4:
                    BigTransScale<int, SepMatrixAccessor<int> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 8:
                    BigTransScale<double, SepMatrixAccessor<double> >( 
                    pMat, rowIndices, colIndices);
                    break;
            }
        }
        else {
            switch (pMat->matrix_type()) {
                case 1:
                    BigTransScale<char, MatrixAccessor<char> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 2:
                    BigTransScale<short, MatrixAccessor<short> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 4:
                    BigTransScale<int, MatrixAccessor<int> >( 
                    pMat, rowIndices, colIndices);
                    break;
                case 8:
                    BigTransScale<double, MatrixAccessor<double> >( 
                    pMat, rowIndices, colIndices);
                    break;
            }
        }
        
        return R_NilValue;
    }
}   // End extern "C"

#endif //BIGEXTENSIONS_UTILS
