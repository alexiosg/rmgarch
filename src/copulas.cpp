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

#include "copulas.h"
using namespace Rcpp;

SEXP copulaNormalC1(SEXP Rbar, SEXP Z)
{
	try {
		Rcpp::NumericMatrix RRbar(Rbar);
		Rcpp::NumericMatrix RZ(Z);
		int n = RZ.nrow(), m = RRbar.nrow(), j;
		List output(2);
		NumericVector Rllh(n);
		arma::mat AR(RRbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		double llhtemp = 0.0, llh = 0.0, temp1 = 0.0, temp2 = 0.0;
		arma::mat AtempR(m, m);
		arma::mat xeye = arma::eye<arma::mat>(m,m);
		AtempR = arma::inv(AR) - xeye;
		temp1 = log(arma::det(AR));
		for(j=0; j<n; j++)
		{
			temp2 = 0.0;
			llhtemp = 0.0;
			temp2 = arma::as_scalar( AZ.row(j) * (AtempR  * arma::trans(AZ.row(j))));
			llhtemp =  0.5 * (temp1 + temp2);
			llh+= llhtemp;
			Rllh[j] = llhtemp;
		}
		output[0] = Rllh;
		output[1] = llh;
	return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->cgarchfit c++ exception (unknown reason)" );
	}
		return R_NilValue;
}

SEXP copulaNormalC2(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP Nbar, SEXP Z, SEXP N, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RN(N);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		List output(4), QtOut(n), RtOut(n);
		NumericVector Rllh(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat AN(RN.begin(), n, m, true);
		double *Rpars = REAL(pars);
		double *Repars = REAL(epars);
		int *Rmodel = INTEGER(model);
		int *Ridx = INTEGER(idx);
		int mo = (int) Repars[2];
		double llhtemp = 0.0, llh = 0.0, temp2 = 0.0;
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat AUQ(m, m);
		arma::mat temp1(m, m);
		arma::mat xeye = arma::eye<arma::mat>(m,m);
		AUQ = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;

		for(j=0; j<mo; j++)
		{
			QtOut[j] = AQbar;
			Rllh[j] = 0.0;
		}
		for(j=mo; j<n; j++)
		{
			temp2 = 0.0;
			llhtemp = 0.0;
			temp1.zeros();
			ARt.zeros();
			AQt = AUQ;
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					AQt = AQt + Rpars[Ridx[1]+i]*( arma::trans( AZ.row(j-(i+1)) ) * AZ.row(j-(i+1)) );
				}
			}
			if(Rmodel[5]>0){
				for(i = 0; i<Rmodel[5]; i++)
				{
					AQt = AQt + Rpars[Ridx[3]+i]*( arma::trans( AN.row(j-(i+1)) ) * AN.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[2]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			arma::vec tmp5 = arma::sqrt(AQt.diag());
			temp1 = tmp5 * arma::trans(tmp5);
			ARt = AQt/temp1;
			RtOut[j] = ARt;
			temp2 = arma::as_scalar( AZ.row(j) * ( (arma::inv(ARt) - xeye) * arma::trans(AZ.row(j))) );
			llhtemp = 0.5*(log(arma::det(ARt)) + temp2);
			llh+= llhtemp;
			Rllh[j] = llhtemp;
		}
		output[0] = QtOut;
		output[1] = Rllh;
		output[2] = llh;
		output[3] = RtOut;
		return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    	} catch(...) {
		::Rf_error( "rmgarch-->cgarchfit c++ exception (unknown reason)" );
    	}
    	return R_NilValue;
}

SEXP copulaStudentC1(SEXP pars, SEXP idx, SEXP Rbar, SEXP Z, SEXP dtZ)
{
	try {
		// dtZ is dt(Z, nu, log = TRUE)
		Rcpp::NumericMatrix RRbar(Rbar);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RdtZ(dtZ);
		int n = RZ.nrow(), m = RRbar.nrow(), j;
		List output(2);
		NumericVector Rllh(n);
		arma::mat AR(RRbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat AdtZ(RdtZ.begin(), n, m, true);
		double *Rpars = REAL(pars);
		int *Ridx = INTEGER(idx);
		double llhtemp = 0.0, llh = 0.0, temp1 = 0.0, temp2 = 0.0, temp4 = 0.0;
		double lcons = Rf_lgammafn(0.5 * (Rpars[Ridx[4]] + m)) - Rf_lgammafn(0.5 * Rpars[Ridx[4]]) - 0.5 * m * log(M_PI*(Rpars[Ridx[4]] - 2));
		arma::mat AtempR(m, m);
		AtempR = arma::inv(AR);
		temp1 = log(arma::det(AR));
		for(j=0; j<n; j++)
		{
			temp2 = 0.0;
			temp4 = 0.0;
			llhtemp = 0.0;
			temp2 = arma::as_scalar( AZ.row(j) * (AtempR * arma::trans(AZ.row(j))));
			temp4 = arma::sum( AdtZ.row(j) );
			llhtemp = lcons - 0.5 * temp1 - 0.5 * (Rpars[Ridx[4]] + m) * log(1.0 + (1.0/(Rpars[Ridx[4]] - 2)) * temp2) - temp4;
			llh+= llhtemp;
			Rllh[j] = -1.0 * llhtemp;
		}
		output[0] = Rllh;
		output[1] = -1.0 * llh;
		return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    	} catch(...) {
		::Rf_error( "rmgarch-->cgarchfit c++ exception (unknown reason)" );
    	}
    	return R_NilValue;
}

SEXP copulaStudentC2(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP Nbar, SEXP Z, SEXP N, SEXP dtZ, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RN(N);
		Rcpp::NumericMatrix RdtZ(dtZ);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		List output(4), QtOut(n), RtOut(n);
		NumericVector Rllh(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat AN(RN.begin(), n, m, true);
		arma::mat AdtZ(RdtZ.begin(), n, m, true);
		// n is already (nrow(data) + maxorder i.e.lag)
		double *Rpars = REAL(pars);
		double *Repars = REAL(epars);
		int *Rmodel = INTEGER(model);
		int *Ridx = INTEGER(idx);
		int mo = (int) Repars[2];
		double llhtemp = 0.0, llh = 0.0, temp2 = 0.0, temp4 = 0.0;
		double lcons = Rf_lgammafn(0.5 * (Rpars[Ridx[4]] + m)) - Rf_lgammafn(0.5 * Rpars[Ridx[4]]) - 0.5 * m * log(M_PI*(Rpars[Ridx[4]] - 2));
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat AUQ(m, m);
		arma::mat temp1(m, m);
		AUQ = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;
		for(j=0; j<mo; j++)
		{
			QtOut[j] = AQbar;
			Rllh[j] = 0.0;
		}
		for(j=mo; j<n; j++)
		{
			temp2 = 0.0;
			temp4 = 0.0;
			llhtemp = 0.0;
			temp1.zeros();
			ARt.zeros();
			AQt = AUQ;
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					AQt = AQt + Rpars[Ridx[1]+i]*( arma::trans( AZ.row(j-(i+1)) ) * AZ.row(j-(i+1)) );
				}
			}
			if(Rmodel[5]>0){
				for(i = 0; i<Rmodel[5]; i++)
				{
					AQt = AQt + Rpars[Ridx[3]+i]*( arma::trans( AN.row(j-(i+1)) ) * AN.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[2]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt/temp1;
			RtOut[j] = ARt;
			temp2 = arma::as_scalar( AZ.row(j) * ( arma::inv(ARt) * arma::trans(AZ.row(j))) );
			temp4 = arma::as_scalar(arma::accu( AdtZ.row(j) ) );
			llhtemp = lcons - 0.5 * log(arma::det(ARt)) - 0.5 * (Rpars[Ridx[4]] + m) * log(1.0 + (1.0/(Rpars[Ridx[4]] - 2)) * temp2) - temp4;
			llh+= llhtemp;
			Rllh[j] = -1.0 * llhtemp;
		}
		output[0] = QtOut;
		output[1] = Rllh;
		output[2] = -1.0 * llh;
		output[3] = RtOut;
		return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    	} catch(...) {
		::Rf_error( "rmgarch-->cgarchfit c++ exception (unknown reason)" );
    	}
    	return R_NilValue;
}

SEXP copuladccsimmvn(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP preQ, SEXP Rbar, SEXP Nbar, SEXP Z, SEXP NZ, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RRbar(Rbar);
		Rcpp::NumericMatrix RNZ(NZ);
		Rcpp::NumericMatrix RpreQ(preQ);
		Rcpp::NumericMatrix RZ(Z);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		List output(3), RtOut(n), QtOut(n), ZtOut(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat ARbar(RRbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat ANZ(RNZ.begin(), n, m, true);
		arma::mat ApreQ(RpreQ.begin(), m, m, true);
		double *Rpars = REAL(pars);
		double *Repars = REAL(epars);
		int *Rmodel = INTEGER(model);
		int *Ridx = INTEGER(idx);
		int mo = (int) Repars[2];
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat temp1(m, m);
		arma::mat AZt(n, m);
		arma::mat ANt(n, m);
		arma::vec eigval(m);
		arma::mat eigvec(m, m);
		arma::mat AUQ(m, m);
		arma::rowvec tmpz(m);
		AUQ = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;
		for(j=0; j<mo; j++)
		{
			QtOut[j] = ApreQ;
			RtOut[j] = ARbar;
			AZt.row(j) = AZ.row(j);
			// Add some small number to avoid division by zero
			ANt.row(j) = ( 0.5*(AZt.row(j)/abs(AZt.row(j)+1E-12)-1) )%AZt.row(j);
		}
		for(j=mo; j<n; j++)
		{
			AQt = AUQ;
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					AQt = AQt + Rpars[Ridx[1]+i]*( arma::trans( AZt.row(j-(i+1)) ) * AZt.row(j-(i+1)) );
				}
			}
			if(Rmodel[5]>0){
				for(i = 0; i<Rmodel[5]; i++)
				{
					AQt = AQt + Rpars[Ridx[3]+i]*( arma::trans( ANt.row(j-(i+1)) ) * ANt.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[2]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			arma::vec tmp5 = arma::sqrt(AQt.diag());
			temp1 = tmp5 * arma::trans(tmp5);
			ARt = AQt/temp1;
			arma::mat Rtmp(m,m);
			// enforce this in case of problems
			RtOut[j] = ARt;
			Rtmp = arma::symmatu(ARt);
			//eig_sym(eigval, eigvec, Rtmp);
			//temp1 = ( eigvec * arma::diagmat( arma::sqrt( eigval ) ) * arma::inv( eigvec ) );
			//AZt.row(j) = ARes.row(j) * temp1;
			// To obtain the negative of the Zs sign(z) = (z/abs(z)-1)/2
			tmpz = ANZ.row(j);
			AZt.row(j) = rmvnormx(Rtmp, tmpz);
			ANt.row(j) = ( 0.5*(AZt.row(j)/abs(AZt.row(j)+1E-12)-1) )%AZt.row(j);
		}
		output[0] = QtOut;
		output[1] = RtOut;
		output[2] = AZt;
		return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->dccsim c++ exception (unknown reason)" );
	}
	return R_NilValue;
}

SEXP copuladccsimmvt(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP preQ, SEXP Rbar, SEXP Nbar, SEXP Z, SEXP NZ, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RRbar(Rbar);
		Rcpp::NumericMatrix RNZ(NZ);
		Rcpp::NumericMatrix RpreQ(preQ);
		Rcpp::NumericMatrix RZ(Z);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		List output(3), RtOut(n), QtOut(n), ZtOut(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat ARbar(RRbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat ANZ(RNZ.begin(), n, m, true);
		arma::mat ApreQ(RpreQ.begin(), m, m, true);
		double *Rpars = REAL(pars);
		double *Repars = REAL(epars);
		int *Rmodel = INTEGER(model);
		int *Ridx = INTEGER(idx);
		int mo = (int) Repars[2];
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat temp1(m, m);
		arma::mat AZt(n, m);
		arma::mat ANt(n, m);
		arma::vec eigval(m);
		arma::mat eigvec(m, m);
		arma::mat AUQ(m, m);
		double nu = Rpars[Ridx[4]];
		arma::rowvec tmpz(m);
		AUQ = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;
		for(j=0; j<mo; j++)
		{
			QtOut[j] = ApreQ;
			RtOut[j] = ARbar;
			AZt.row(j) = AZ.row(j);
			// Add some small number to avoid division by zero
			ANt.row(j) = ( 0.5*(AZt.row(j)/abs(AZt.row(j)+1E-12)-1) )%AZt.row(j);
		}
		for(j=mo; j<n; j++)
		{
			AQt = AUQ;
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					AQt = AQt + Rpars[Ridx[1]+i]*( arma::trans( AZt.row(j-(i+1)) ) * AZt.row(j-(i+1)) );
				}
			}
			if(Rmodel[5]>0){
				for(i = 0; i<Rmodel[5]; i++)
				{
					AQt = AQt + Rpars[Ridx[3]+i]*( arma::trans( ANt.row(j-(i+1)) ) * ANt.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[2]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			arma::vec tmp5 = arma::sqrt(AQt.diag());
			temp1 = tmp5 * arma::trans(tmp5);
			ARt = AQt/temp1;
			arma::mat Rtmp(m,m);
			// enforce this in case of problems
			RtOut[j] = ARt;
			Rtmp = arma::symmatu(ARt);
			tmpz = ANZ.row(j);
			AZt.row(j) = rmvtx(Rtmp, nu, tmpz);
			// To obtain the negative of the Zs sign(z) = (z/abs(z)-1)/2
			ANt.row(j) = ( 0.5*(AZt.row(j)/abs(AZt.row(j)+1E-12)-1) )%AZt.row(j);
		}
		output[0] = QtOut;
		output[1] = RtOut;
		output[2] = AZt;
		return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "rmgarch-->dccsim c++ exception (unknown reason)" );
	}
	return R_NilValue;
}
