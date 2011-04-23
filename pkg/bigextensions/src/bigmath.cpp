#ifndef BIGEXTENSIONS_MATH
#define BIGEXTENSIONS_MATH

#include "bigextensions/bigex.h"
#include <math.h>

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
void BigLog(BigMatrix *pMat) {
    BIG_FUNCTION_SETUP(BMAccessorType)
    BIG_FUNCTION_LOOP_EVAL(double, CType, log)
    return;
}

template<typename CType, typename BMAccessorType>
void BigLog10(BigMatrix *pMat) {
    BIG_FUNCTION_SETUP(BMAccessorType)
    BIG_FUNCTION_LOOP_EVAL(double, CType, log10)
    return;
}

template<typename CType, typename BMAccessorType>
void BigAbs(BigMatrix *pMat) {
    BIG_FUNCTION_SETUP(BMAccessorType)
    BIG_FUNCTION_LOOP_EVAL(double, CType, fabs)
    return;
}

template<typename CType, typename BMAccessorType>
void BigSqrt(BigMatrix *pMat) {
    BIG_FUNCTION_SETUP(BMAccessorType)
    BIG_FUNCTION_LOOP_EVAL(double, CType, sqrt)
    return;
}

template<typename CType, typename BMAccessorType>
void BigPow(BigMatrix *pMat, SEXP exponent) {
    BIG_FUNCTION_SETUP(BMAccessorType)
    
    double value = NUMERIC_DATA(exponent)[0];
    
    CType *pColumn;    \
    for (i = 0; i < numCols; ++i) {   \
        pColumn = m[i];   \
        for (j = 0; j < numRows; ++j)     \
            if (!ISNAN(pColumn[j]))   \
                pColumn[j] = (CType)(pow((double)pColumn[j], value));  \
    }
    
    return;
}


extern "C" {
    
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
    

    SEXP BigTanhMain(SEXP bigaddr) {       
        CALL_CPLUSPLUS_FUN_NOVALUE(BigTanh)
        return R_NilValue;
    }

    SEXP BigAtanhMain(SEXP bigaddr) {
        CALL_CPLUSPLUS_FUN_NOVALUE(BigAtanh)
        return R_NilValue;
    }
    
    SEXP BigLogMain(SEXP bigaddr) {
        CALL_CPLUSPLUS_FUN_NOVALUE(BigLog)
        return R_NilValue;
    }
    
    SEXP BigLog10Main(SEXP bigaddr) {
        CALL_CPLUSPLUS_FUN_NOVALUE(BigLog10)
        return R_NilValue;
    }
    
    SEXP BigAbsMain(SEXP bigaddr) {
        CALL_CPLUSPLUS_FUN_NOVALUE(BigAbs)
        return R_NilValue;
    }
    
    SEXP BigSqrtMain(SEXP bigaddr) {
        CALL_CPLUSPLUS_FUN_NOVALUE(BigSqrt)
        return R_NilValue;
    }
    
    
    #define CALL_CPLUSPLUS_FUN_ONEVALUE(fun) \
      BigMatrix *pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(bigaddr)); \
      if (pMat->separated_columns()) {  \
          switch (pMat->matrix_type()) {    \
              case 1:   \
                  fun<char, SepMatrixAccessor<char> >(   \
                  pMat, arg);     \
                  break;    \
              case 2:   \
                  fun<short, SepMatrixAccessor<short> >(     \
                  pMat, arg);    \
                  break;    \
              case 4:   \
                  fun<int, SepMatrixAccessor<int> >(     \
                  pMat, arg);  \
                  break;    \
              case 8:   \
                  fun<double, SepMatrixAccessor<double> >(   \
                  pMat, arg);     \
                  break;    \
          }     \
      }     \
      else {    \
          switch (pMat->matrix_type()) {    \
              case 1:   \
                  fun<char, MatrixAccessor<char> >(      \
                  pMat, arg);     \
                  break;    \
              case 2:   \
                  fun<short, MatrixAccessor<short> >(    \
                  pMat, arg);    \
                  break;    \
              case 4:   \
                  fun<int, MatrixAccessor<int> >(    \
                  pMat, arg);  \
                  break;    \
              case 8:   \
                  fun<double, MatrixAccessor<double> >(      \
                  pMat, arg);     \
                  break;    \
          }     \
      }
    
      SEXP BigPowMain(SEXP bigaddr, SEXP arg) {
          CALL_CPLUSPLUS_FUN_ONEVALUE(BigPow)
          return R_NilValue;
      }

}

#endif //BIGEXTENSIONS_MATH
