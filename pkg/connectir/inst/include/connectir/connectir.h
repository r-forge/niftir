#ifndef _connectir_CONNECTIR
#define _connectir_CONNECTIR

#include <Rcpp.h>
#include <RcppArmadillo.h>
//using namespace Rcpp; 
//using namespace Rcpp::sugar;

#include "connectirDefines.h"

#include <math.h>
#include <iostream>

#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/bigmemoryDefines.h"
#include "bigmemory/isna.hpp"

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

extern "C" {
    SEXP ComputePvalsMain(SEXP Rinmat, SEXP Routmat, SEXP Routcol);
}

RcppExport SEXP CombineSubMapsMain(SEXP LIST_allVoxs_allSubs, SEXP ADDR_oneVox_allSubs, SEXP Rseed_index, SEXP Rvoxindices, SEXP Rnvoxs, SEXP Rnsubs);

RcppExport SEXP CombineSubMapsTransSimpleMain(SEXP LIST_allVoxs_allSubs, SEXP ADDR_oneVox_allSubs, SEXP Rseed_index, SEXP Rvoxindices, SEXP Rnvoxs, SEXP Rnsubs);

// utils.cpp
BigMatrix* sbm_to_bm(SEXP Sbm);
void sbm_to_arma_xd(SEXP SbM, arma::mat& M);
void sub_sbm_to_arma_xd(SEXP SbM, arma::mat& M, SEXP SfirstCol, SEXP SlastCol);
RcppExport SEXP test_func(SEXP Sbm);

// arith.cpp
RcppExport SEXP big_add_multiply_scalar(SEXP SX, SEXP SX,  
                                        SEXP Sa, SEXP Sb, 
                                        SEXP SX_firstCol, SEXP SX_lastCol, 
                                        SEXP SY_firstCol, SEXP SY_lastCol);

// summary_stats.cpp
RcppExport SEXP big_rowsum(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol);
RcppExport SEXP big_rowmean(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol);
RcppExport SEXP big_colsum(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol);
RcppExport SEXP big_colmean(SEXP SX, SEXP SX_firstCol, SEXP SX_lastCol);

// qlm.cpp
RcppExport SEXP big_qlm_rank(SEXP Xr);
RcppExport SEXP big_qlm_dd(SEXP Xr);
RcppExport SEXP big_qlm_fit(SEXP yr, SEXP Xr, SEXP coefr, SEXP residr, SEXP mser, SEXP nr, SEXP kr, SEXP mr);
RcppExport SEXP big_qlm_residuals(SEXP yr, SEXP Xr, SEXP residr);
RcppExport SEXP big_qlm_contrasts(SEXP fit_coefr, SEXP fit_mser, SEXP conr, SEXP ddr, SEXP coefr, SEXP ser, SEXP tvalr, SEXP mr);

// cor.cpp
RcppExport SEXP test_sub_matrix(SEXP As, SEXP As_firstCol, SEXP As_lastCol);
RcppExport SEXP big_cor(SEXP As, SEXP Bs, SEXP Cs, SEXP A_firstCol, SEXP A_lastCol, SEXP B_firstCol, SEXP B_lastCol, SEXP C_firstCol, SEXP C_lastCol);
RcppExport SEXP big_tcor(SEXP As, SEXP Bs, SEXP Cs, SEXP A_firstCol, SEXP A_lastCol, SEXP B_firstCol, SEXP B_lastCol, SEXP C_firstCol, SEXP C_lastCol);

// subdist2.cpp
RcppExport SEXP big_gower(SEXP SX, SEXP SY,  
                          SEXP SX_firstCol, SEXP SX_lastCol, 
                          SEXP SY_firstCol, SEXP SY_lastCol);

#endif // _connectir_CONNECTIR
