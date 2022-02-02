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
#include "rmdist.h"

arma::rowvec rmvnormx(arma::mat R, arma::rowvec Z){
	Rcpp::RNGScope scope;
	int m = R.n_rows;
	arma::vec eigval(m);
	arma::mat eigvec(m, m);
	arma::mat temp(m, m);
	arma::eig_sym(eigval, eigvec, R);
	arma::rowvec ans(m);
	temp = ( eigvec * arma::diagmat( arma::sqrt( eigval ) ) * arma::inv( eigvec ) );
	ans = Z * temp;
	return ans;
}

arma::rowvec rmvtx(arma::mat R, const double nu, arma::rowvec Z){
	Rcpp::RNGScope scope;
	int m = R.n_rows;
	double rc = Rf_rchisq(nu);
	arma::mat RR = ((nu-2.0)/nu)*R;
	double v = sqrt(nu/rc);
	arma::rowvec ans(m);
	ans = v*rmvnormx(RR, Z);
	return ans;
}

arma::rowvec rmvlx(arma::mat R, arma::rowvec Z){
	Rcpp::RNGScope scope;
	int m = R.n_rows;
	double e  = Rf_rexp(1);
	arma::rowvec ans(m);
	ans = sqrt(e) * rmvnormx(R, Z);
	return ans;
}
