#include "connectir/connectir.h"

// In R:
// library(connectir)
// x <- as.big.matrix(matrix(1:20, 5, 4), type="double")
// y <- .Call("test_sub_matrix", x, as.double(2), as.double(4))
// all.equal(x[,2:4], y)
SEXP test_sub_matrix(SEXP As, SEXP As_firstCol, SEXP As_lastCol) {
    try {
        BM_COL_INIT(As_firstCol, A_firstCol)
        BM_COL_INIT(As_lastCol, A_lastCol)
    
        BM_TO_ARMA_INIT()
        
        SUB_BM_TO_ARMA_MULTIPLE(As, A, A_firstCol, A_lastCol)
        
        return Rcpp::wrap( A );
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

// Cs <- (t(As) %*% Bs)/(nrow(Bs)-1)
SEXP big_cor(SEXP As, SEXP Bs, SEXP Cs, 
             SEXP As_firstCol, SEXP As_lastCol, 
             SEXP Bs_firstCol, SEXP Bs_lastCol, 
             SEXP Cs_firstCol, SEXP Cs_lastCol) 
{    
    try {
        BM_COL_INIT(As_firstCol, A_firstCol)
        BM_COL_INIT(As_lastCol, A_lastCol)
        BM_COL_INIT(Bs_firstCol, B_firstCol)
        BM_COL_INIT(Bs_lastCol, B_lastCol)
        BM_COL_INIT(Cs_firstCol, C_firstCol)
        BM_COL_INIT(Cs_lastCol, C_lastCol)
        
        BM_TO_ARMA_INIT()
        SUB_BM_TO_ARMA_MULTIPLE(As, A, A_firstCol, A_lastCol)
        SUB_BM_TO_ARMA_MULTIPLE(Bs, B, B_firstCol, B_lastCol)
        double df = pMat->nrow() - 1;
        SUB_BM_TO_ARMA_MULTIPLE(Cs, C, C_firstCol, C_lastCol)
        
        // TODO: make this an explicit user option!
        if (ncol == 1 && B.n_cols != 1 && C.n_rows == (A.n_cols*B.n_cols))
            C.reshape(A.n_cols, B.n_cols);
        
        if (A.n_cols != C.n_rows)
            Rf_error("incorrect number of rows in output matrix");
        if (B.n_cols != C.n_cols)
            Rf_error("incorrect number of columns in output matrix");
        
        C = (arma::trans(A) * B)/df;
        
        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}

// Cs <- (As %*% t(Bs))/(nrow(As)-1)
SEXP big_tcor(SEXP As, SEXP Bs, SEXP Cs, 
             SEXP As_firstCol, SEXP As_lastCol, 
             SEXP Bs_firstCol, SEXP Bs_lastCol, 
             SEXP Cs_firstCol, SEXP Cs_lastCol) 
{    
    try {
        BM_COL_INIT(As_firstCol, A_firstCol)
        BM_COL_INIT(As_lastCol, A_lastCol)
        BM_COL_INIT(Bs_firstCol, B_firstCol)
        BM_COL_INIT(Bs_lastCol, B_lastCol)
        BM_COL_INIT(Cs_firstCol, C_firstCol)
        BM_COL_INIT(Cs_lastCol, C_lastCol)
        
        BM_TO_ARMA_INIT()
        SUB_BM_TO_ARMA_MULTIPLE(As, A, A_firstCol, A_lastCol)
        SUB_BM_TO_ARMA_MULTIPLE(Bs, B, B_firstCol, B_lastCol)
        double df = ncol - 1;
        SUB_BM_TO_ARMA_MULTIPLE(Cs, C, C_firstCol, C_lastCol)
        
        // TODO: make this an explicit user option!
        if (ncol == 1 && B.n_rows != 1 && C.n_rows == (A.n_rows*B.n_rows))
            C.reshape(A.n_rows, B.n_rows);
        
        if (A.n_rows != C.n_rows)
            Rf_error("incorrect number of rows in output matrix");
        if (B.n_rows != C.n_cols)
            Rf_error("incorrect number of columns in output matrix");
        
        C = (A * arma::trans(B))/df;
        
        return R_NilValue;
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    
    return R_NilValue;
}

double kendall_worker_d(SEXP SRatings)
{
    try {
        arma::mat Ratings(1,1);
        const double* old_rptr = sbm_to_arma_xd(SRatings, Ratings);
    
        arma::uvec Ovec(Ratings.n_rows);
        arma::mat Rmat(Ratings.n_rows, Ratings.n_cols);
        double* Rvec;
    
        // Ranking
        arma::u32 i, j;
        for (i = 0; i < Ratings.n_cols; ++i) {
            Ovec = arma::sort_index(Ratings.unsafe_col(i));
            Rvec = Rmat.colptr(i);
            for (j = 0; j < Ratings.n_rows; ++j) {
                Rvec[Ovec(j)] = static_cast<double>(j+1);
            }
        }
    
        // rowSums
        arma::colvec rs = arma::sum(Rmat, 1);
    
        // var
        double v = var(rs);
    
        // clean up
        Rvec = NULL;
        free_arma(Ratings, old_rptr);
    
        return v;
    } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return -1;
}

SEXP kendall_worker(SEXP SRatings)
{
    try {
        return Rcpp::wrap( kendall_worker_d(SRatings) );
    } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue;
}

SEXP voxelwise_kendall(SEXP Slist_CorMaps, SEXP SSeedMaps, SEXP Snseeds)
{
    try {
        Rcpp::List list_CorMaps(Slist_CorMaps);
        index_type nsubs = list_CorMaps.size();
        index_type nseeds = static_cast<index_type>(DOUBLE_DATA(Snseeds)[0]);
        
        SEXP SSubMaps;
        arma::mat SubMaps(1,1);
        const double* old_sptr = SubMaps.memptr();
        
        arma::mat SeedMaps(1,1);
        const double* old_smptr = sbm_to_arma_xd(SSeedMaps, SeedMaps);
        index_type nvoxs = static_cast<index_type>(SeedMaps.n_rows);
        if (nsubs != (index_type)SeedMaps.n_cols)
            Rf_error("number of columns in seedmap incorrect");
        
        arma::vec ks(nseeds);
        
        index_type seedi, subi, voxi;
        double* smap;
        // Loop through each seed voxel
        for (seedi = 0; seedi < nseeds; ++seedi) {
            
            // Combine seed maps across subjects for given voxel
            for (subi = 0; subi < nsubs; ++subi)
            {
                PROTECT(SSubMaps = VECTOR_ELT(Slist_CorMaps, subi));
                sbm_to_arma_xd(SSubMaps, SubMaps);
                UNPROTECT(1);
                smap = SeedMaps.colptr(subi);
                for (voxi = 0; voxi < nvoxs; ++voxi) {
                    smap[voxi] = SubMaps(seedi,voxi);
                }
            }
            
            // Get kendall's w
            ks(seedi) = kendall_worker_d(SSeedMaps);
        }
        
        // Scale
        double d_nvoxs = static_cast<double>(nvoxs);
        double d_nsubs = static_cast<double>(nsubs);
        ks = (12*ks*d_nvoxs-1) / (pow(d_nsubs,2)*(pow(d_nvoxs,3)-d_nvoxs));
        
        smap = NULL;
        free_arma(SubMaps, old_sptr);
        free_arma(SeedMaps, old_smptr);
        
        return Rcpp::wrap( ks );
    } catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) { 
      ::Rf_error( "c++ exception (unknown reason)" ); 
    }
    return R_NilValue;
}
