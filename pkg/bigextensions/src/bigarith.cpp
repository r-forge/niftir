#ifndef BIGEXTENSIONS_ARITH
#define BIGEXTENSIONS_ARITH

#include "bigextensions/bigex.h"


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
              pColumn[j] = (TYPE)(value OPERATOR (double)pColumn[j]);  \
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


extern "C" {
    
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
    
    
    #define CALL_CPLUSPLUS_COLS_FUN(fun) \
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
      
      // Do some operation between each value in a column and another value

      SEXP BigAddColsMain(SEXP bigaddr, SEXP values, SEXP narm) {       
          CALL_CPLUSPLUS_COLS_FUN(BigAddCols)
          return R_NilValue;
      }

      SEXP BigSubtractColsMain(SEXP bigaddr, SEXP values, SEXP narm) {
          CALL_CPLUSPLUS_COLS_FUN(BigSubtractCols)
          return R_NilValue;
      }

      SEXP BigMultiplyColsMain(SEXP bigaddr, SEXP values, SEXP narm) {
          CALL_CPLUSPLUS_COLS_FUN(BigMultiplyCols)
          return R_NilValue;
      }

      SEXP BigDivideColsMain(SEXP bigaddr, SEXP values, SEXP narm) {
          CALL_CPLUSPLUS_COLS_FUN(BigDivideCols)
          return R_NilValue;
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
      
      // Do some operation between one scalar value and every value in the matrix

      SEXP BigAddScalarMain(SEXP bigaddr, SEXP values, SEXP narm) {       
          CALL_CPLUSPLUS_FUN(BigAddScalar)
          return R_NilValue;
      }

      SEXP BigSubtractScalarMain(SEXP bigaddr, SEXP values, SEXP narm) {
          CALL_CPLUSPLUS_FUN(BigSubtractScalar)
          return R_NilValue;
      }
}

#endif //BIGEXTENSIONS_ARITH
