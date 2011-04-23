#ifndef BIGEXTENSIONS_ANALYTICS_PLUS
#define BIGEXTENSIONS_ANALYTICS_PLUS

#include "bigextensions/bigex.h"


template<typename CType, typename BMAccessorType>
void BigSum(BigMatrix *mat, SEXP colIndices, SEXP rowIndices, double *value) {
    BMAccessorType m( *mat );
    
    double *cols = NUMERIC_DATA(colIndices);
    double *rows = NUMERIC_DATA(rowIndices);
    
    index_type i=0;
    index_type j=0;
    index_type xj=0;
    index_type numCols = GET_LENGTH(colIndices);
    index_type numRows = GET_LENGTH(rowIndices);
    CType *pColumn;
    
    double s = 0;
    
    for (i = 0; i < numCols; ++i) {
        pColumn = m[static_cast<index_type>(cols[i])-1];
        for (j = 0; j < numRows; ++j) {
            xj = static_cast<index_type>(rows[j])-1;
            s += (double)pColumn[xj];
        }
    }
    
    value[0] = s;
    
    return;
}

#define BIG_FUNCTION_COLAPPLY_SETUP(ACCESSOR_TYPE, COL_TYPE)        \
  ACCESSOR_TYPE m( *pMat );                                         \
                                                                    \
  LDOUBLE s = 0.0;                                                  \
  index_type naCount=0;                                             \
                                                                    \
  Rboolean Rnarm = (Rboolean)LOGICAL_VALUE(narm);                   \
  index_type i=0;                                                   \
  index_type j=0;                                                   \
  index_type curj=0;                                                \
  COL_TYPE *pColumn;


template<typename CType, typename BMAccessorType>
void BigColMean(BigMatrix *pMat, double *pCols, index_type nCols, 
                double *pRows, index_type nRows, double *pRet, SEXP narm) 
{
    BIG_FUNCTION_COLAPPLY_SETUP(BMAccessorType, CType);
    
    for (i=0; i < nCols; ++i) {
        s = 0.0;
        naCount = 0;
        pColumn = m[static_cast<index_type>(pCols[i])-1];
        
        for (j=0; j < nRows; ++j) {
            curj = static_cast<index_type>(pRows[j])-1;
            if (!ISNAN(pColumn[curj])) {
                s += (LDOUBLE)pColumn[curj];
            } else if (!Rnarm) {
                naCount = nRows;
                break;
            } else {
                ++naCount;
            }
        }
        
        if (nRows-naCount > 0)
            s /= (LDOUBLE)(nRows-naCount);
        else
            s = NA_REAL;
        
        pRet[i] = (double)s;
    }
    
    return;
}

template<typename CType, typename BMAccessorType>
void BigColZMean(BigMatrix *pMat, double *pCols, index_type nCols, 
                double *pRows, index_type nRows, double *pRet, SEXP narm) 
{
    BIG_FUNCTION_COLAPPLY_SETUP(BMAccessorType, CType);
    
    for (i=0; i < nCols; ++i) {
        s = 0.0;
        naCount = 0;
        pColumn = m[static_cast<index_type>(pCols[i])-1];
        
        for (j=0; j < nRows; ++j) {
            curj = static_cast<index_type>(pRows[j])-1;
            if (!ISNAN(pColumn[curj])) {
                s += (LDOUBLE)(atanh(pColumn[curj]));
            } else if (!Rnarm) {
                naCount = nRows;
                break;
            } else {
                ++naCount;
            }
        }
        
        if (nRows-naCount > 0)
            s /= (LDOUBLE)(nRows-naCount);
        else
            s = NA_REAL;
        
        pRet[i] = (LDOUBLE)(tanh(s));
    }
    
    return;
}


extern "C" {

    SEXP BigSumMain(SEXP addr, SEXP cols, SEXP rows) {
        SEXP ret = R_NilValue;
        ret = PROTECT(NEW_NUMERIC(1));
        double *pRet = NUMERIC_DATA(ret);
        
        BigMatrix *pMat = (BigMatrix*)R_ExternalPtrAddr(addr);
        
        if (pMat->separated_columns()) {
            switch (pMat->matrix_type()) {
                case 1:
                    BigSum<char, SepMatrixAccessor<char> >(
                    pMat, cols, rows, pRet);
                    break;
                case 2:
                    BigSum<short, SepMatrixAccessor<short> >(
                    pMat, cols, rows, pRet);
                    break;
                case 4:
                    BigSum<int, SepMatrixAccessor<int> >(
                    pMat, cols, rows, pRet);
                    break;
                case 8:
                    BigSum<double, SepMatrixAccessor<double> >(
                    pMat, cols, rows, pRet);
                    break;
            }
        }
        else {
            switch (pMat->matrix_type()) {
                case 1:
                    BigSum<char, MatrixAccessor<char> >(
                    pMat, cols, rows, pRet);
                    break;
                case 2:
                    BigSum<short, MatrixAccessor<short> >(
                    pMat, cols, rows, pRet);
                    break;
                case 4:
                    BigSum<int, MatrixAccessor<int> >(
                    pMat, cols, rows, pRet);
                    break;
                case 8:
                    BigSum<double, MatrixAccessor<double> >(
                    pMat, cols, rows, pRet);
                    break;
            }
        }
        
        UNPROTECT(1);
        return(ret);
    }
    
    #define colapply_setup(fun)                                             \
      BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(addr));\
      double *pCols = NUMERIC_DATA(cols);                                   \
      index_type nCols = GET_LENGTH(cols);                                  \
      double *pRows = NUMERIC_DATA(rows);                                   \
      index_type nRows = GET_LENGTH(rows);                                  \
                                                                            \
      SEXP ret = R_NilValue;                                                \
      ret = PROTECT(NEW_NUMERIC(nCols));                                    \
      double *pRet = NUMERIC_DATA(ret);                                     \
                                                                            \
      if (pMat->separated_columns()) {                                      \
          switch (pMat->matrix_type()) {                                    \
              case 1:                                                       \
                  fun<char, SepMatrixAccessor<char> >(                      \
                  pMat, pCols, nCols, pRows, nRows, pRet, narm);            \
                  break;                                                    \
              case 2:                                                       \
                  fun<short, SepMatrixAccessor<short> >(                    \
                  pMat, pCols, nCols, pRows, nRows, pRet, narm);            \
                  break;                                                    \
              case 4:                                                       \
                  fun<int, SepMatrixAccessor<int> >(                        \
                  pMat, pCols, nCols, pRows, nRows, pRet, narm);            \
                  break;                                                    \
              case 8:                                                       \
                  fun<double, SepMatrixAccessor<double> >(                  \
                  pMat, pCols, nCols, pRows, nRows, pRet, narm);            \
                  break;                                                    \
          }                                                                 \
      }                                                                     \
      else {                                                                \
          switch (pMat->matrix_type()) {                                    \
              case 1:                                                       \
                  fun<char, MatrixAccessor<char> >(                         \
                  pMat, pCols, nCols, pRows, nRows, pRet, narm);            \
                  break;                                                    \
              case 2:                                                       \
                  fun<short, MatrixAccessor<short> >(                       \
                  pMat, pCols, nCols, pRows, nRows, pRet, narm);            \
                  break;                                                    \
              case 4:                                                       \
                  fun<int, MatrixAccessor<int> >(                           \
                  pMat, pCols, nCols, pRows, nRows, pRet, narm);            \
                  break;                                                    \
              case 8:                                                       \
                  fun<double, MatrixAccessor<double> >(                     \
                  pMat, pCols, nCols, pRows, nRows, pRet, narm);            \
                  break;                                                    \
          }                                                                 \
      }

    SEXP BigColMeanMain(SEXP addr, SEXP cols, SEXP rows, SEXP narm) {
      colapply_setup(BigColMean);
      UNPROTECT(1);
      return(ret);
    }

    SEXP BigColZMeanMain(SEXP addr, SEXP cols, SEXP rows, SEXP narm) {
      colapply_setup(BigColZMean);
      UNPROTECT(1);
      return(ret);
    }
}

#endif //BIGEXTENSIONS_ANALYTICS_PLUS
