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

#include "ica.h"
using namespace Rcpp;

SEXP radicalrot(SEXP X, SEXP idx)
{
	try {
		Rcpp::NumericMatrix RX(X);
		int *ridx = INTEGER(idx);
		int m = ridx[0];
		int k = ridx[1];
		int nx = RX.nrow(), mx = RX.ncol(), i;
		arma::mat AX(RX.begin(), nx, mx, true);
		arma::mat rot(2,2);
		arma::mat rotp(nx, mx);
		double theta = 0.0;
		arma::rowvec ent(k);
		arma::rowvec vals(mx);
		arma::rowvec margtheta(2);
		arma::rowvec tmp;
		for(i=0;i<k;i++){
			theta = ( (double) i/((double) k-1)*(0.5*M_PI)-(0.25*M_PI));
			rot.zeros();
			rotp.zeros();
			rot(0,0) = cos(theta);
			rot(0,1) = -1.0*sin(theta);
			rot(1,0) = sin(theta);
			rot(1,1) = cos(theta);
			rotp = rot * AX;
			vals = rotp.row(0);
			vals = arma::sort(vals, 0);
			tmp = arma::log(vals.cols(m,mx-1) - vals.cols(0,mx-m-1));
			margtheta(0) = arma::as_scalar(arma::accu(tmp));
			vals = rotp.row(1);
			vals = arma::sort(vals, 0);
			tmp = arma::log(vals.cols(m,mx-1) - vals.cols(0,mx-m-1));
			margtheta(1) = arma::as_scalar(arma::accu(tmp));
			ent(i) = arma::as_scalar(arma::accu(margtheta));
		}
		return wrap(ent);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "africa-->radical c++ exception (unknown reason)" );
	}
		return R_NilValue;
}

arma::mat fs10(arma::mat X, arma::mat B, const double nsamples)
{
	//Rcpp::NumericMatrix RX(X);
	//Rcpp::NumericMatrix RB(B);
	//int nrx = RX.nrow(), ncx = RX.ncol();
	//int nrb = RB.nrow(), ncb = RB.ncol();
	//arma::mat AX(RX.begin(), nrx, ncx, true);
	//arma::mat AB(RB.begin(), nrb, ncb, true);
	//double *rv = REAL(v);
	//double nsamples = rv[0];
	arma::mat A = (X * arma::pow(X.t()*B, 3))/nsamples - 3*B;
	return(A);
}

arma::mat fs11(arma::mat X, arma::mat B, const double nsamples, const double myy)
{
	arma::mat Y  = X.t()*B;
	arma::mat Gpow3 = arma::pow(Y, 3);
	arma::mat Beta =  arma::sum(Y%Gpow3);
	arma::mat D = 1/(Beta -3*nsamples);
	arma::mat A = B + myy*B*(Y.t()*Gpow3 - arma::diagmat(Beta)) * arma::diagmat(D);
	return(A);
}

arma::mat sqrtm(arma::mat X)
{
	arma::mat U;
	arma::vec s;
	arma::mat V;
	arma::svd(U,s,V,X);
	arma::mat A = U * arma::sqrt(arma::diagmat(s)) * U.t();
	return(A);
}
