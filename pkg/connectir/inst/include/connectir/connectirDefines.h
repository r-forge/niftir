#ifndef _connectir_DEFINES_H
#define _connectir_DEFINES_H

// take an 'input' R big matrix variable to an 'output' double armadillo matrix
#define BM_TO_ARMA(INPUT, OUTPUT) \
    Rcpp::RObject BM_TO_ARMA_bm(INPUT); \
    SEXP BM_TO_ARMA_addr = BM_TO_ARMA_bm.slot("address"); \
    BigMatrix *BM_TO_ARMA_pMat = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(BM_TO_ARMA_addr)); \
    if (BM_TO_ARMA_pMat->matrix_type() != 8) \
        Rf_error("Big Matrix must be of type double"); \
    double *BM_TO_ARMA_ptr_double = reinterpret_cast<double*>(BM_TO_ARMA_pMat->matrix()); \
    arma::mat OUTPUT(BM_TO_ARMA_ptr_double, BM_TO_ARMA_pMat->nrow(), BM_TO_ARMA_pMat->ncol(), false);

#define SET_ACCESSOR(ptr, mat)                                          \
    if (ptr->separated_columns()) {                                     \
        switch (ptr->matrix_type()) {                                   \
            case 1:                                                     \
                SepMatrixAccessor<char> mat( *ptr );                    \
                break;                                                  \
            case 2:                                                     \
                SepMatrixAccessor<short> mat( *ptr );                   \
                break;                                                  \
            case 4:                                                     \
                SepMatrixAccessor<int> mat( *ptr );                     \
                break;                                                  \
            case 8:                                                     \
                SepMatrixAccessor<double> mat( *ptr );                  \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (ptr->matrix_type()) {                                   \
            case 1:                                                     \
                MatrixAccessor<char> mat( *ptr );                       \
                break;                                                  \
            case 2:                                                     \
                MatrixAccessor<short> mat( *ptr );                      \
                break;                                                  \
            case 4:                                                     \
                MatrixAccessor<int> mat( *ptr );                        \
                break;                                                  \
            case 8:                                                     \
                MatrixAccessor<double> mat( *ptr );                     \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_ONE(FUN, INBIGMAT)                        \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT);                                              \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT);                                              \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT);                                              \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT);                                              \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT);                                              \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT);                                              \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT);                                              \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT);                                              \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_TWO(FUN, INBIGMAT, ARG2)                  \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT, ARG2);                                        \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_THREE(FUN, INBIGMAT, ARG2, ARG3)                  \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT, ARG2, ARG3);                                        \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_FOUR(FUN, INBIGMAT, ARG2, ARG3, ARG4)                  \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT, ARG2, ARG3, ARG4);                                        \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_FIVE(FUN, INBIGMAT, ARG2, ARG3, ARG4, ARG5)                  \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5);                                        \
                break;                                                  \
        }                                                               \
    }

#define CALL_BIGFUNCTION_ARGS_SIX(FUN, INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6)                  \
    SEXP ret;                                                           \
    if (INBIGMAT->separated_columns()) {                                \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, SepMatrixAccessor<char> >(                    \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, SepMatrixAccessor<short> >(                  \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, SepMatrixAccessor<int> >(                      \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, SepMatrixAccessor<double> >(                \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
        }                                                               \
    }                                                                   \
    else {                                                              \
        switch (INBIGMAT->matrix_type()) {                              \
            case 1:                                                     \
                ret = FUN<char, MatrixAccessor<char> >(                       \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 2:                                                     \
                ret = FUN<short, MatrixAccessor<short> >(                     \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 4:                                                     \
                ret = FUN<int, MatrixAccessor<int> >(                         \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
            case 8:                                                     \
                ret = FUN<double, MatrixAccessor<double> >(                   \
                INBIGMAT, ARG2, ARG3, ARG4, ARG5, ARG6);                                        \
                break;                                                  \
        }                                                               \
    }


#endif // _connectir_DEFINES_H
