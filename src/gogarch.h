/*################################################################################
##
##   R package rmgarch by Alexios Galanos Copyright (C) 2009-2013
##   This file is part of the R package rmgarch.
##
##   The R package rmgarch is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package rmgarch is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
#################################################################################*/
#ifndef _gogarch_H
#define _gogarch_H
#include <RcppArmadillo.h>
RcppExport SEXP tvbetacovar(SEXP wi, SEXP Vi, SEXP di);
RcppExport SEXP tvbetacoskew(SEXP wi, SEXP Si, SEXP di);
RcppExport SEXP tvbetacokurt(SEXP wi, SEXP Ki, SEXP di);
RcppExport SEXP gogarchSigma(SEXP S, SEXP A);
RcppExport SEXP gogarchCov(SEXP S, SEXP A);
RcppExport SEXP gogarchCor(SEXP S, SEXP A);
RcppExport SEXP Cov2Cor(SEXP , SEXP );
RcppExport SEXP ArrayDiag(SEXP , SEXP );
RcppExport SEXP Cov2Res(SEXP , SEXP , SEXP );
RcppExport SEXP RowApplySort(SEXP );
RcppExport SEXP RowUnique(SEXP );
RcppExport SEXP gogarchcssigma(SEXP , SEXP );
RcppExport SEXP gogarchcksigma(SEXP , SEXP );
arma::mat coskewind( arma::rowvec );
#endif
