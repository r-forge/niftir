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
RcppExport SEXP bm_rowsum(SEXP Rbigmat);
RcppExport SEXP bm_rowmean(SEXP Rbigmat);

#endif // _connectir_CONNECTIR
