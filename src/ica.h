/*################################################################################
##
##   R package rmgarch by Alexios Ghalanos Copyright (C) 2009-2013
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
#ifndef _ica_H
#define _ica_H
#include <RcppArmadillo.h>
RcppExport SEXP radicalrot(SEXP , SEXP);
arma::mat fs10(arma::mat X, arma::mat B, const double nsamples);
arma::mat fs11(arma::mat X, arma::mat B, const double nsamples, const double myy);
arma::mat sqrtm(arma::mat X);
#endif
