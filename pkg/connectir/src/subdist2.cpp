#include "connectir/connectir.h"

//// 
//SEXP big_subdist_worker(SEXP Slist_cormaps, SEXP Sdmat, INDS?, VOXS?)
//{
//    // 1. Create matrix of subject connectivity maps at one node/voxel
//    
//    // 2. Copy over subject connectivity maps from list to matrix
//    //    and scale
//    
//    // 3. Correlate
//    
//    // 4. Convert to distance
//}
//
//.Call("CombineSubMapsMain", sub.cormaps, subsMap@address, as.double(i), 
//      as.double(voxs[-cor_inds[i]]), as.double(nvoxs-1), as.double(nsubs))
//big_cor(x=subsMap, z=outmat, z_firstCol=dist_inds[i], z_lastCol=dist_inds[i])
//.Call("big_add_scalar", outmat, as.double(-1), as.double(1), 
//        as.double(dist_inds[i]), as.double(dist_inds[i]));


// Y <- -0.5 * dmat^2 %*% (I - ones %*% t(ones)/n)
SEXP big_gower(SEXP SX, SEXP SY, 
               SEXP SX_firstCol, SEXP SX_lastCol, 
               SEXP SY_firstCol, SEXP SY_lastCol)
{
    try {
        arma::mat X;
        sub_sbm_to_arma_xd(SX, X, SX_firstCol, SX_lastCol);
        arma::mat Y;
        sub_sbm_to_arma_xd(SY, Y, SY_firstCol, SY_lastCol);
                
        if (X.n_rows != Y.n_rows || X.n_cols != Y.n_cols)
            Rf_error("dimension mismatch between input and output matrices");
        double n = sqrt(static_cast<double>(Y.n_rows));
        
        using namespace arma;
        
        // adj = I - ones %*% t(ones)/n
        mat I = eye<mat>(n, n);
        colvec z_ones = ones<colvec>(n);
        mat adj = I - z_ones * trans(z_ones) / n;
        
        // A = -0.5 * dmat^2
        Y = -(arma::pow(X, 2)/2);
        
        // G = A %*% adj
        mat Y_vox;
        for (index_type i = 0; i < Y.n_rows; ++i)
        {
            Y_vox = mat(Y.colptr(i), n, n, false);  // This actually creates a copy
            Y_vox = Y_vox * adj;
            Y_vox.set_size(Y.n_rows, 1);
            Y.col(i) = Y_vox;
        }
        
        return Rcpp::wrap( SY ); // or return R_NilValue?
    } catch(std::exception &ex) {
        forward_exception_to_r(ex);
    } catch(...) {
        ::Rf_error("c++ exception (unknown reason)");
    }
    return R_NilValue;
}
