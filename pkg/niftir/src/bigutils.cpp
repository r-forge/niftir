#include <iostream> // hack to make sure we are using the right "length" 
                    // function
#include <algorithm>

#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/bigmemoryDefines.h"
#include "bigmemory/isna.hpp"

#include <math.h>
#include <R.h>
#include <Rdefines.h>

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

#define BIG_MAT_SETUP(ACCESSOR_TYPE)    \
  ACCESSOR_TYPE xm( *pXmat );     \
  ACCESSOR_TYPE ym( *pYmat );     \
  ACCESSOR_TYPE zm( *pZmat );     \
    \
  double *xCols = NUMERIC_DATA(xColIndices);   \
  double *yCols = NUMERIC_DATA(yColIndices);   \
  double *zCols = NUMERIC_DATA(zColIndices);   \
  double *xRows = NUMERIC_DATA(xRowIndices);   \
  double *yRows = NUMERIC_DATA(yRowIndices);   \
  double *zRows = NUMERIC_DATA(zRowIndices);   \
    \
  index_type i=0;   \
  index_type j=0;   \
  index_type numCols = GET_LENGTH(zColIndices);    \
  index_type numRows = GET_LENGTH(zRowIndices);   \
  index_type xc=0;  \
  index_type yc=0;  \
  index_type zc=0;  \
  index_type xr=0;  \
  index_type yr=0;  \
  index_type zr=0;

#define BIG_MAT_LOOP_EVAL(TYPE, OPERATOR) \
  for (i = 0; i < numCols; ++i) {   \
      xc = static_cast<index_type>(xCols[i])-1;  \
      yc = static_cast<index_type>(yCols[i])-1;  \
      zc = static_cast<index_type>(zCols[i])-1;  \
      for (j = 0; j < numRows; ++j) {     \
          xr = static_cast<index_type>(xRows[j])-1;  \
          yr = static_cast<index_type>(yRows[j])-1;  \
          zr = static_cast<index_type>(zRows[j])-1;  \
          if (ISNAN(xm[xc][xr]) || ISNAN(ym[yc][yr]))   \
              zm[zc][zr] = NA_VALUE;  \
          else \
              zm[zc][zr] = (TYPE)xm[xc][xr] OPERATOR (TYPE)ym[yc][yr];   \
      } \
  }

template<typename CType, typename BMAccessorType>
void BigAddMats(BigMatrix *pXmat, BigMatrix *pYmat, BigMatrix *pZmat, SEXP xColIndices, SEXP yColIndices, SEXP zColIndices, SEXP xRowIndices, SEXP yRowIndices, SEXP zRowIndices, CType NA_VALUE) {
    BIG_MAT_SETUP(BMAccessorType)
    BIG_MAT_LOOP_EVAL(CType, +)
    return;
}

template<typename CType, typename BMAccessorType>
void BigSubtractMats(BigMatrix *pXmat, BigMatrix *pYmat, BigMatrix *pZmat, SEXP xColIndices, SEXP yColIndices, SEXP zColIndices, SEXP xRowIndices, SEXP yRowIndices, SEXP zRowIndices, CType NA_VALUE) {
    BIG_MAT_SETUP(BMAccessorType)   
    BIG_MAT_LOOP_EVAL(CType, -)
    return;
}

template<typename CType, typename BMAccessorType>
void BigMultiplyMats(BigMatrix *pXmat, BigMatrix *pYmat, BigMatrix *pZmat, SEXP xColIndices, SEXP yColIndices, SEXP zColIndices, SEXP xRowIndices, SEXP yRowIndices, SEXP zRowIndices, CType NA_VALUE) {
    BIG_MAT_SETUP(BMAccessorType)   
    BIG_MAT_LOOP_EVAL(CType, *)
    return;
}

template<typename CType, typename BMAccessorType>
void BigDivideMats(BigMatrix *pXmat, BigMatrix *pYmat, BigMatrix *pZmat, SEXP xColIndices, SEXP yColIndices, SEXP zColIndices, SEXP xRowIndices, SEXP yRowIndices, SEXP zRowIndices, CType NA_VALUE) {
    BIG_MAT_SETUP(BMAccessorType)
    BIG_MAT_LOOP_EVAL(CType, /)
    return;
}


#define BIG_COL_SETUP(ACCESSOR_TYPE)    \
  ACCESSOR_TYPE m( *pMat );     \
  Rboolean narm = (Rboolean)LOGICAL_VALUE(Rnarm);   \
    \
  index_type i=0;   \
  index_type j=0;   \
  index_type numCols = pMat->ncol();    \
  index_type numRows = pMat->nrow();    \
    \
  double *pColValues = NUMERIC_DATA(Rvalues);   \
  index_type numValues = GET_LENGTH(Rvalues);   \
  if (numCols != numValues)     \
    error("number of values does not match number of columns");     \

#define BIG_COL_LOOP_EVAL(TYPE, OPERATOR) \
  TYPE *pColumn;    \
  double curvalue;  \
  for (i = 0; i < numCols; ++i) {   \
      pColumn = m[i];   \
      curvalue = pColValues[i];     \
      for (j = 0; j < numRows; ++j)     \
          if (!ISNAN(pColumn[j]))   \
              pColumn[j] = (TYPE)((double)pColumn[j] OPERATOR curvalue);  \
  }

template<typename CType, typename BMAccessorType>
void BigAddCols(BigMatrix *pMat, SEXP Rvalues, SEXP Rnarm, CType NA_VALUE) {
    BIG_COL_SETUP(BMAccessorType)
    BIG_COL_LOOP_EVAL(CType, +)
    return;
}

template<typename CType, typename BMAccessorType>
void BigSubtractCols(BigMatrix *pMat, SEXP Rvalues, SEXP Rnarm, CType NA_VALUE) {
    BIG_COL_SETUP(BMAccessorType)   
    BIG_COL_LOOP_EVAL(CType, -)
    return;
}

template<typename CType, typename BMAccessorType>
void BigMultiplyCols(BigMatrix *pMat, SEXP Rvalues, SEXP Rnarm, CType NA_VALUE) {
    BIG_COL_SETUP(BMAccessorType)   
    BIG_COL_LOOP_EVAL(CType, *)
    return;
}

template<typename CType, typename BMAccessorType>
void BigDivideCols(BigMatrix *pMat, SEXP Rvalues, SEXP Rnarm, CType NA_VALUE) {
    BIG_COL_SETUP(BMAccessorType)
    BIG_COL_LOOP_EVAL(CType, /)
    return;
}

  
#define BIG_SCALAR_SETUP(ACCESSOR_TYPE)    \
  ACCESSOR_TYPE m( *pMat );     \
  Rboolean narm = (Rboolean)LOGICAL_VALUE(Rnarm);   \
    \
  double value = NUMERIC_DATA(Rvalues)[0];   \
    \
  index_type i=0;   \
  index_type j=0;   \
  index_type numCols = pMat->ncol();    \
  index_type numRows = pMat->nrow();


#define BIG_SCALAR_LOOP_EVAL(TYPE, OPERATOR) \
  TYPE *pColumn;    \
  for (i = 0; i < numCols; ++i) {   \
      pColumn = m[i];   \
      for (j = 0; j < numRows; ++j)     \
          if (!ISNAN(pColumn[j]))   \
              pColumn[j] = (TYPE)((double)pColumn[j] OPERATOR value);  \
  }

template<typename CType, typename BMAccessorType>
void BigAddScalar(BigMatrix *pMat, SEXP Rvalues, SEXP Rnarm, CType NA_VALUE) {
    BIG_SCALAR_SETUP(BMAccessorType)
    BIG_SCALAR_LOOP_EVAL(CType, +)
    return;
}

template<typename CType, typename BMAccessorType>
void BigSubtractScalar(BigMatrix *pMat, SEXP Rvalues, SEXP Rnarm, CType NA_VALUE) {
    BIG_SCALAR_SETUP(BMAccessorType)   
    BIG_SCALAR_LOOP_EVAL(CType, -)
    return;
}

template<typename CType, typename BMAccessorType>
void BigMultiplyScalar(BigMatrix *pMat, SEXP Rvalues, SEXP Rnarm, CType NA_VALUE) {
    BIG_SCALAR_SETUP(BMAccessorType)   
    BIG_SCALAR_LOOP_EVAL(CType, *)
    return;
}

template<typename CType, typename BMAccessorType>
void BigDivideScalar(BigMatrix *pMat, SEXP Rvalues, SEXP Rnarm, CType NA_VALUE) {
    BIG_SCALAR_SETUP(BMAccessorType)
    BIG_SCALAR_LOOP_EVAL(CType, /)
    return;
}

template<typename CType, typename BMAccessorType>
void BigPow(BigMatrix *pMat) {
    BMAccessorType m( *pMat );
    
    index_type i=0;
    index_type j=0;
    index_type numCols = pMat->ncol();
    index_type numRows = pMat->nrow();
    
    CType *pColumn;    \
    for (i = 0; i < numCols; ++i) {   \
        pColumn = m[i];   \
        for (j = 0; j < numRows; ++j)     \
            if (!ISNAN(pColumn[j]))   \
                pColumn[j] = (CType)(pow((double)pColumn[j], 2));  \
    }
    
    return;
}

#define BIG_FUNCTION_SETUP(ACCESSOR_TYPE)    \
  ACCESSOR_TYPE m( *pMat );     \
    \
  index_type i=0;   \
  index_type j=0;   \
  index_type numCols = pMat->ncol();    \
  index_type numRows = pMat->nrow();


#define BIG_FUNCTION_LOOP_EVAL(IN_TYPE, OUT_TYPE, FUNCTION) \
  OUT_TYPE *pColumn;    \
  for (i = 0; i < numCols; ++i) {   \
      pColumn = m[i];   \
      for (j = 0; j < numRows; ++j)     \
          if (!ISNAN(pColumn[j]))   \
              pColumn[j] = (OUT_TYPE)(FUNCTION((IN_TYPE)pColumn[j]));  \
  }

template<typename CType, typename BMAccessorType>
void BigTanh(BigMatrix *pMat) {
    BIG_FUNCTION_SETUP(BMAccessorType)
    BIG_FUNCTION_LOOP_EVAL(double, CType, tanh)
    return;
}

template<typename CType, typename BMAccessorType>
void BigAtanh(BigMatrix *pMat) {
    BIG_FUNCTION_SETUP(BMAccessorType)
    BIG_FUNCTION_LOOP_EVAL(double, CType, atanh)
    return;
}

template<typename CType, typename BMAccessorType>
void BigSqrt(BigMatrix *pMat) {
    BIG_FUNCTION_SETUP(BMAccessorType)
    BIG_FUNCTION_LOOP_EVAL(double, CType, sqrt)
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
    
//    index_type rowOffset = pOrigMat->row_offset();
//    index_type colOffset = pOrigMat->col_offset();
    
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
    
   #define CALL_CPLUSPLUS_MATS_FUN(fun) \
     BigMatrix *pXmat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(xaddr)); \
     BigMatrix *pYmat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(yaddr)); \
     BigMatrix *pZmat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(zaddr)); \
     if (pXmat->separated_columns() != pYmat->separated_columns() && \
         pXmat->separated_columns() != pZmat->separated_columns()) \
           error("all big matrices are not the same column separated type"); \
     if (pZmat->separated_columns()) {  \
         switch (pZmat->matrix_type()) {    \
             case 1:   \
                 fun<char, SepMatrixAccessor<char> >(   \
                 pXmat, pYmat, pZmat, xCols, yCols, zCols, xRows, yRows, zRows, \
                 NA_CHAR);     \
                 break;    \
             case 2:   \
                 fun<short, SepMatrixAccessor<short> >(     \
                 pXmat, pYmat, pZmat, xCols, yCols, zCols, xRows, yRows, zRows, \
                 NA_SHORT);     \
                 break;    \
             case 4:   \
                 fun<int, SepMatrixAccessor<int> >(     \
                 pXmat, pYmat, pZmat, xCols, yCols, zCols, xRows, yRows, zRows, \
                 NA_INTEGER);     \
                 break;    \
             case 8:   \
                 fun<double, SepMatrixAccessor<double> >(   \
                 pXmat, pYmat, pZmat, xCols, yCols, zCols, xRows, yRows, zRows, \
                 NA_REAL);     \
                 break;    \
         }     \
     }     \
     else {    \
         switch (pZmat->matrix_type()) {    \
             case 1:   \
                 fun<char, MatrixAccessor<char> >(      \
                 pXmat, pYmat, pZmat, xCols, yCols, zCols, xRows, yRows, zRows, \
                 NA_CHAR);     \
                 break;    \
             case 2:   \
                 fun<short, MatrixAccessor<short> >(    \
                 pXmat, pYmat, pZmat, xCols, yCols, zCols, xRows, yRows, zRows, \
                 NA_SHORT);     \
                 break;    \
             case 4:   \
                 fun<int, MatrixAccessor<int> >(    \
                 pXmat, pYmat, pZmat, xCols, yCols, zCols, xRows, yRows, zRows, \
                 NA_INTEGER);     \
                 break;    \
             case 8:   \
                 fun<double, MatrixAccessor<double> >(      \
                 pXmat, pYmat, pZmat, xCols, yCols, zCols, xRows, yRows, zRows, \
                 NA_REAL);     \
                 break;    \
         }     \
     }
    
#define CALL_CPLUSPLUS_FUN(fun) \
  BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(bigaddr)); \
  if (pMat->separated_columns()) {  \
      switch (pMat->matrix_type()) {    \
          case 1:   \
              fun<char, SepMatrixAccessor<char> >(   \
              pMat, values, narm, NA_CHAR);     \
              break;    \
          case 2:   \
              fun<short, SepMatrixAccessor<short> >(     \
              pMat, values, narm, NA_SHORT);    \
              break;    \
          case 4:   \
              fun<int, SepMatrixAccessor<int> >(     \
              pMat, values, narm, NA_INTEGER);  \
              break;    \
          case 8:   \
              fun<double, SepMatrixAccessor<double> >(   \
              pMat, values, narm, NA_REAL);     \
              break;    \
      }     \
  }     \
  else {    \
      switch (pMat->matrix_type()) {    \
          case 1:   \
              fun<char, MatrixAccessor<char> >(      \
              pMat, values, narm, NA_CHAR);     \
              break;    \
          case 2:   \
              fun<short, MatrixAccessor<short> >(    \
              pMat, values, narm, NA_SHORT);    \
              break;    \
          case 4:   \
              fun<int, MatrixAccessor<int> >(    \
              pMat, values, narm, NA_INTEGER);  \
              break;    \
          case 8:   \
              fun<double, MatrixAccessor<double> >(      \
              pMat, values, narm, NA_REAL);     \
              break;    \
      }     \
  }
  
 #define CALL_CPLUSPLUS_FUN_NOVALUE(fun) \
   BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(bigaddr)); \
   if (pMat->separated_columns()) {  \
       switch (pMat->matrix_type()) {    \
           case 1:   \
               fun<char, SepMatrixAccessor<char> >(   \
               pMat);     \
               break;    \
           case 2:   \
               fun<short, SepMatrixAccessor<short> >(     \
               pMat);    \
               break;    \
           case 4:   \
               fun<int, SepMatrixAccessor<int> >(     \
               pMat);  \
               break;    \
           case 8:   \
               fun<double, SepMatrixAccessor<double> >(   \
               pMat);     \
               break;    \
       }     \
   }     \
   else {    \
       switch (pMat->matrix_type()) {    \
           case 1:   \
               fun<char, MatrixAccessor<char> >(      \
               pMat);     \
               break;    \
           case 2:   \
               fun<short, MatrixAccessor<short> >(    \
               pMat);    \
               break;    \
           case 4:   \
               fun<int, MatrixAccessor<int> >(    \
               pMat);  \
               break;    \
           case 8:   \
               fun<double, MatrixAccessor<double> >(      \
               pMat);     \
               break;    \
       }     \
   }
    
 
    // Do some operation between the same ith,jth element of two big matrices 
    
    SEXP BigAddMatsMain(SEXP xaddr, SEXP yaddr, SEXP zaddr, SEXP xCols, SEXP yCols, SEXP zCols, SEXP xRows, SEXP yRows, SEXP zRows, SEXP narm) {       
        CALL_CPLUSPLUS_MATS_FUN(BigAddMats)
        return R_NilValue;
    }
    
    SEXP BigSubtractMatsMain(SEXP xaddr, SEXP yaddr, SEXP zaddr, SEXP xCols, SEXP yCols, SEXP zCols, SEXP xRows, SEXP yRows, SEXP zRows, SEXP narm) {
        CALL_CPLUSPLUS_MATS_FUN(BigSubtractMats)
        return R_NilValue;
    }
    
    SEXP BigMultiplyMatsMain(SEXP xaddr, SEXP yaddr, SEXP zaddr, SEXP xCols, SEXP yCols, SEXP zCols, SEXP xRows, SEXP yRows, SEXP zRows, SEXP narm) {
        CALL_CPLUSPLUS_MATS_FUN(BigMultiplyMats)
        return R_NilValue;
    }
    
    SEXP BigDivideMatsMain(SEXP xaddr, SEXP yaddr, SEXP zaddr, SEXP xCols, SEXP yCols, SEXP zCols, SEXP xRows, SEXP yRows, SEXP zRows, SEXP narm) {
        CALL_CPLUSPLUS_MATS_FUN(BigDivideMats)
        return R_NilValue;
    }
    
    
    // Do some operation between each value in a column and another value
    
    SEXP BigAddColsMain(SEXP bigaddr, SEXP values, SEXP narm) {       
        CALL_CPLUSPLUS_FUN(BigAddCols)
        return R_NilValue;
    }
    
    SEXP BigSubtractColsMain(SEXP bigaddr, SEXP values, SEXP narm) {
        CALL_CPLUSPLUS_FUN(BigSubtractCols)
        return R_NilValue;
    }
    
    SEXP BigMultiplyColsMain(SEXP bigaddr, SEXP values, SEXP narm) {
        CALL_CPLUSPLUS_FUN(BigMultiplyCols)
        return R_NilValue;
    }
    
    SEXP BigDivideColsMain(SEXP bigaddr, SEXP values, SEXP narm) {
        CALL_CPLUSPLUS_FUN(BigDivideCols)
        return R_NilValue;
    }
    
    // Do some operation between one scalar value and every value in the matrix
    
    SEXP BigAddScalarMain(SEXP bigaddr, SEXP values, SEXP narm) {       
        CALL_CPLUSPLUS_FUN(BigAddScalar)
        return R_NilValue;
    }
    
    SEXP BigSubtractScalarMain(SEXP bigaddr, SEXP values, SEXP narm) {
        CALL_CPLUSPLUS_FUN(BigSubtractScalar)
        return R_NilValue;
    }
    
    SEXP BigMultiplyScalarMain(SEXP bigaddr, SEXP values, SEXP narm) {
        CALL_CPLUSPLUS_FUN(BigMultiplyScalar)
        return R_NilValue;
    }
    
    SEXP BigDivideScalarMain(SEXP bigaddr, SEXP values, SEXP narm) {
        CALL_CPLUSPLUS_FUN(BigDivideScalar)
        return R_NilValue;
    }
    
    // Apply some function to each value in the matrix
    
    SEXP BigTanhMain(SEXP bigaddr) {       
        CALL_CPLUSPLUS_FUN_NOVALUE(BigTanh)
        return R_NilValue;
    }
    
    SEXP BigAtanhMain(SEXP bigaddr) {
        CALL_CPLUSPLUS_FUN_NOVALUE(BigAtanh)
        return R_NilValue;
    }
    
    SEXP BigPowMain(SEXP bigaddr) {
        CALL_CPLUSPLUS_FUN_NOVALUE(BigPow)
        return R_NilValue;
    }
    
    SEXP BigSqrtMain(SEXP bigaddr) {
        CALL_CPLUSPLUS_FUN_NOVALUE(BigSqrt)
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
    
    
    
    // Mean...
    
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
    
}   // End extern "C"
