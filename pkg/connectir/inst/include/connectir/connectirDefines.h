#ifndef _connectir_DEFINES_H
#define _connectir_DEFINES_H

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
