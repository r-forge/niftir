//#ifndef _connectir_SUBDIST
//#define _connectir_SUBDIST

#include "connectirDefines.h"

#include <Rcpp.h>
using namespace Rcpp; 
using namespace Rcpp::sugar;

#include <math.h>
#include <iostream>

#include "bigmemory/BigMatrix.h"
#include "bigmemory/MatrixAccessor.hpp"
#include "bigmemory/bigmemoryDefines.h"
#include "bigmemory/isna.hpp"

//#include <R.h>
//#include <Rinternals.h>
//#include <Rdefines.h>

RcppExport SEXP CombineSubMapsMain(SEXP LIST_allVoxs_allSubs, SEXP ADDR_oneVox_allSubs, SEXP Rseed_index, SEXP Rvoxindices, SEXP Rnvoxs, SEXP Rnsubs);

//#endif // _connectir_SUBDIST
