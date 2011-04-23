#include "connectir/subdist.h"

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

#define Z_STOP \
    BEGIN_RCPP \
    throw std::range_error("debugging error"); \
    END_RCPP


template<typename CType, typename BMAccessorType>
SEXP CombineSubMaps(BigMatrix *oneVox_allSubs, SEXP RallVoxs_allSubs, index_type seed, double *pVoxs, index_type nvoxs, index_type nsubs) {
    
    List allVoxs_allSubs(RallVoxs_allSubs);
    BMAccessorType outMat( *oneVox_allSubs );
    
    //BEGIN_RCPP
    if (nsubs != oneVox_allSubs->ncol())
        Rf_error("nsubs must equal oneVox_allSubs->ncol");
    if (nvoxs != oneVox_allSubs->nrow())
        Rf_error("nsubs must equal oneVox_allSubs->ncol");
    //END_RCPP
    
    // loop through each subject's map
    index_type s = 0;
    index_type v = 0;
    index_type vv = 0;
    LDOUBLE x = 0;
    LDOUBLE delta = 0;
    LDOUBLE mean = 0;
    LDOUBLE M2 = 0;
    LDOUBLE stdev = 0;
    SEXP Rp;
//    CType *inCol;
    CType *outCol;
    LDOUBLE scaled_x;
    BigMatrix *allVoxs_oneSub;
    SEXP tmp;
    for (s=0; s < nsubs; ++s) {
        tmp = allVoxs_allSubs[s];
        RObject RallVoxs_oneSub(tmp);
        Rp = RallVoxs_oneSub.slot("address");
        allVoxs_oneSub = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(Rp));
        BMAccessorType inMat( *allVoxs_oneSub );
        
//        inCol = inMat[seed];
        delta = mean = M2 = stdev = 0;
        for (v=0; v < nvoxs; ++v) {
            // todo: add checking for NaN...but shouldn't really have any!
            // maybe can also pass the exact list of voxs to loop through!
            // if (!ISNAN(pColumn[curj]))
            // NA_REAL
            vv = static_cast<index_type>(pVoxs[v]-1);
            x = static_cast<LDOUBLE>(inMat[vv][seed]);
            delta = x - mean;
            mean = mean + delta/static_cast<LDOUBLE>(v+1);
            M2 = M2 + delta*(x - mean);
        }
        stdev = sqrt(M2/(static_cast<LDOUBLE>(nvoxs-1)));
        //printf("mean: %f; stdev: %f\n", mean, stdev);
        
        outCol = outMat[s];
        for (v=0; v < nvoxs; ++v) {
            vv = static_cast<index_type>(pVoxs[v]-1);
            scaled_x = (static_cast<LDOUBLE>(inMat[vv][seed])-mean)/stdev;
            outCol[v] = static_cast<CType>(scaled_x);
        }
    }
    
    return R_NilValue;
}


SEXP CombineSubMapsMain(SEXP LIST_allVoxs_allSubs, SEXP ADDR_oneVox_allSubs, SEXP Rseed_index, SEXP Rvoxindices, SEXP Rnvoxs, SEXP Rnsubs) {
        BigMatrix *oneVox_allSubs = reinterpret_cast<BigMatrix*>(R_ExternalPtrAddr(ADDR_oneVox_allSubs));
        index_type seed = static_cast<index_type>(NUMERIC_DATA(Rseed_index)[0]) - 1;
        double *pVoxs = NUMERIC_DATA(Rvoxindices);
        index_type nvoxs = static_cast<index_type>(NUMERIC_DATA(Rnvoxs)[0]);
        index_type nsubs = static_cast<index_type>(NUMERIC_DATA(Rnsubs)[0]);
        
        CALL_BIGFUNCTION_ARGS_SIX(CombineSubMaps, oneVox_allSubs, LIST_allVoxs_allSubs, seed, pVoxs, nvoxs, nsubs)
        
        return(ret);
}

