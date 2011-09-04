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
RcppExport SEXP big_add_scalar_right(SEXP Sbigmat, SEXP Svalue, SEXP Sicol, SEXP Sncol);
RcppExport SEXP big_add_scalar_left(SEXP Sbigmat, SEXP Svalue, SEXP Sicol, SEXP Sncol);
RcppExport SEXP bm_rowsum(SEXP Rbigmat);
RcppExport SEXP bm_rowmean(SEXP Rbigmat);

// qlm.cpp
RcppExport SEXP big_qlm_rank(SEXP Xr);
RcppExport SEXP big_qlm_dd(SEXP Xr);
RcppExport SEXP big_qlm_fit(SEXP yr, SEXP Xr, SEXP coefr, SEXP residr, SEXP mser, SEXP nr, SEXP kr, SEXP mr);
RcppExport SEXP big_qlm_residuals(SEXP yr, SEXP Xr, SEXP residr);
RcppExport SEXP big_qlm_contrasts(SEXP fit_coefr, SEXP fit_mser, SEXP conr, SEXP ddr, SEXP coefr, SEXP ser, SEXP tvalr, SEXP mr);

// cor.cpp
RcppExport SEXP big_cor(SEXP As, SEXP Bs, SEXP Cs, SEXP Cicol, SEXP Cncol);
RcppExport SEXP big_tcor(SEXP As, SEXP Bs, SEXP Cs, SEXP Cicol, SEXP Cncol);
RcppExport SEXP big_ztransform(SEXP Cs);

#endif // _connectir_CONNECTIR
