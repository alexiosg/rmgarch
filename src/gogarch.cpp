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
// Extra method for speedup
#include "gogarch.h"
#undef CS
using namespace Rcpp;

SEXP tvbetacovar(SEXP wi, SEXP Vi, SEXP di)
{
	try {
	NumericVector CV(Vi);
	int *d = INTEGER(di);
	arma::cube V(CV.begin(), d[0], d[1], d[2]);
	arma::mat w = as<arma::mat>(wi);
	arma::vec A(d[2]);
	int i;
	for(i=0;i<d[2];i++){
		arma::mat temp = V.slice(i);
		A(i) = arma::as_scalar(w.row(i) * arma::trans(temp.row(d[3])))/arma::as_scalar(w.row(i) * temp * arma::trans(w.row(i))) ;
	}
	return wrap( A );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->GOGARCH beta covariance extractor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP tvbetacoskew(SEXP wi, SEXP Si, SEXP di)
{
	try {
	NumericVector CS(Si);
	int *d = INTEGER(di);
	arma::cube S(CS.begin(), d[0], d[1], d[2]);
	arma::mat w = as<arma::mat>(wi);
	arma::vec A(d[2]);
	int i;
	for(i=0;i<d[2];i++){
		arma::vec ww = arma::trans(arma::kron(w.row(i), w.row(i)));
		arma::mat stemp = S.slice(i);
		A(i) = arma::as_scalar(stemp.row(d[3]) * ww)/arma::as_scalar(w.row(i) *stemp * ww);
	}
	return wrap( A );
	} catch( std::exception &ex ) {
	forward_exception_to_r( ex );
	} catch(...) {
	::Rf_error( "rmgarch-->GOGARCH beta coskew extractor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP tvbetacokurt(SEXP wi, SEXP Ki, SEXP di)
{
	try {
	NumericVector CK(Ki);
	int *d = INTEGER(di);
	arma::cube K(CK.begin(), d[0], d[1], d[2]);
	arma::mat w = as<arma::mat>(wi);
	arma::vec A(d[2]);
	int i;
	for(i=0;i<d[2];i++){
		arma::vec ww = arma::trans(arma::kron(w.row(i), arma::kron(w.row(i), w.row(i))));
		arma::mat ktemp = K.slice(i);
		A(i) = arma::as_scalar(ktemp.row(d[3]) * ww)/arma::as_scalar(w.row(i) *ktemp * ww);
	}
	return wrap( A );
	} catch( std::exception &ex ) {
	forward_exception_to_r( ex );
	} catch(...) {
	::Rf_error( "rmgarch-->GOGARCH beta cokurt extractor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP gogarchSigma(SEXP S, SEXP A)
{
	try{
		Rcpp::NumericMatrix RS(S);
		Rcpp::NumericMatrix RA(A);
		int m = RS.ncol(), n = RS.nrow(), mm = RA.ncol(), mn = RA.nrow(), i;
		arma::mat AS(RS.begin(), n, m, true);
		arma::mat AA(RA.begin(), mn, mm, true);
		arma::mat S(n, mn);
		for(i=0;i<n;i++){
			S.row(i) = arma::trans(arma::diagvec(AA * arma::diagmat(AS.row(i))*AA.t()));
		}
		return wrap( S );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->gogarch extractor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP gogarchCov(SEXP S, SEXP A)
{
	try{
		Rcpp::NumericMatrix RS(S);
		Rcpp::NumericMatrix RA(A);
		int m = RS.ncol(), n = RS.nrow(), mm = RA.ncol(), mn = RA.nrow(), i;
		arma::mat AS(RS.begin(), n, m, true);
		arma::mat AA(RA.begin(), mn, mm, true);
		arma::cube S(mn, mn, n);
		for(i=0;i<n;i++){
			S.slice(i) = AA * arma::diagmat(AS.row(i))*AA.t();
		}
		return wrap( S );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->gogarch extractor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP gogarchCor(SEXP S, SEXP A)
{
	try{
		Rcpp::NumericMatrix RS(S);
		Rcpp::NumericMatrix RA(A);
		int m = RS.ncol(), n = RS.nrow(), mm = RA.ncol(), mn = RA.nrow(), i;
		arma::mat AS(RS.begin(), n, m, true);
		arma::mat AA(RA.begin(), mn, mm, true);
		arma::cube S(mn, mn, n);
		for(i=0;i<n;i++){
			arma::mat tmp1 = AA * arma::diagmat(AS.row(i))*AA.t();
			arma::mat tmp2 = arma::diagmat(1/arma::sqrt(diagvec(tmp1)));
			S.slice(i) = tmp2 * tmp1 * tmp2.t();
		}
		return wrap( S );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->gogarch extractor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}


SEXP Cov2Cor(SEXP YY, SEXP dimm){
	try{
		Rcpp::NumericMatrix Y(YY);
		Rcpp::IntegerVector dim(dimm);
		arma::cube AY(Y.begin(), dim[0], dim[1], dim[2]);
		arma::cube AX(dim[0], dim[1], dim[2]);
		int i;
		arma::mat tmp(dim[0], dim[1]);
		for(i=0;i<dim[2];i++){
			tmp = arma::diagmat(1/arma::sqrt(diagvec(AY.slice(i))));
			AX.slice(i) = tmp * AY.slice(i) * tmp.t();
		}
		return wrap( AX );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->Cov2Cor c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP ArrayDiag(SEXP YY, SEXP dimm){
	try{
		int i;
		Rcpp::NumericVector Y(YY);
		Rcpp::IntegerVector dim(dimm);
		arma::cube AY(Y.begin(), dim[0], dim[1], dim[2]);
		arma::mat AX(dim[1],dim[2]);
		for(i=0;i<dim[2];i++){
			AX.col(i) =  arma::diagvec(AY.slice(i));
		}
		AX = AX.t();
		return(wrap(AX));
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
	::Rf_error( "rmgarch-->ArrayDiag c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP Cov2Res(SEXP YY, SEXP ZZ, SEXP dimm){
	try{
		int i;
		Rcpp::NumericVector Y(YY);
		Rcpp::IntegerVector dim(dimm);
		Rcpp::NumericMatrix Z(ZZ);
		arma::cube AY(Y.begin(), dim[1], dim[1], dim[0]);
		arma::mat AZ(Z.begin(), dim[0], dim[1]);
		arma::mat AR(dim[0], dim[1]);
		for(i=0;i<dim[0];i++){
			AR.row(i) =  AZ.row(i) * diagmat(arma::sqrt(diagvec(AY.slice(i))));
		}
		return(wrap(AR));
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
	::Rf_error( "rmgarch-->Cov2Res c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP RowApplySort(SEXP YY){
	try{
		Rcpp::NumericMatrix Y(YY);
		int m = Y.ncol(), n = Y.nrow(), i;
		arma::mat AY(Y.begin(), n, m, true);
		arma::mat AX(n, m);
		for(i=0;i<n;i++){
			AX.row(i) = sort(AY.row(i));
		}
		return wrap( AX );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->RowApplySort c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP RowUnique(SEXP YY){
	try{
		Rcpp::NumericMatrix Y(YY);
		int m = Y.ncol(), n = Y.nrow(), i;
		arma::mat AY(Y.begin(), n, m, true);
		arma::mat AX(n, 2);
		for(i=0;i<n;i++){
			AX.row(i) = unique(AY.row(i));
		}
		return wrap( AX );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->RowUnique c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP gogarchcssigma(SEXP idx, SEXP SS)
{
	try{
		Rcpp::NumericMatrix ix(idx);
		int m = ix.nrow(), i;
		arma::mat yix(ix.begin(), m, 3, false);
		arma::vec S = Rcpp::as<arma::vec>(SS);
		arma::vec R(m);
		arma::uvec dx(3);
		for(i=0;i<m;i++){
			dx(0) = (int) yix(i,0);
			dx(1) = (int) yix(i,1);
			dx(2) = (int) yix(i,2);
			R(i) = arma::prod(S.elem(dx));
		}
		return wrap( R );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->gogarchcssigma c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP gogarchcksigma(SEXP idx, SEXP SS)
{
	try{
		Rcpp::NumericMatrix ix(idx);
		int m = ix.nrow(), i;
		arma::mat yix(ix.begin(), m, 4, false);
		arma::vec S = Rcpp::as<arma::vec>(SS);
		arma::vec R(m);
		arma::uvec dx(4);
		for(i=0;i<m;i++){
			dx(0) = (int) yix(i,0);
			dx(1) = (int) yix(i,1);
			dx(2) = (int) yix(i,2);
			dx(3) = (int) yix(i,3);
			R(i) = arma::prod(S.elem(dx));
		}
		return wrap( R );
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->gogarchcksigma c++ exception (unknown reason)" );
	}
	return R_NilValue;
}


arma::mat coskewind(arma::rowvec skew){
	int n = skew.n_cols, i;
	int nn = n*n;
	arma::mat S(n, nn);
	S.zeros();
	for(i=0;i<n;i++){
		int indx = (i)*nn + (i)*n + i;
		S(indx) = skew(i);
	}
	return(S);
}
