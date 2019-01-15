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
#ifndef _rmdist_H
#define _rmdist_H
#undef trunc
#undef ftrunc
#define R_NO_REMAP
#define STRICT_R_HEADERS
#include <R.h>
#include <RcppArmadillo.h>
arma::rowvec rmvnormx(arma::mat , arma::rowvec );
arma::rowvec rmvtx(arma::mat , const double , arma::rowvec );
arma::rowvec rmvlx(arma::mat , arma::rowvec );
#endif
