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

#include "rdcc.h"
using namespace Rcpp;
/* Rcpp inputs appended with R*
 * ARMA inputs appended with A*
 */
SEXP dccnormC1(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP Nbar, SEXP Z, SEXP N, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RN(N);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		List output(3), QtOut(n);
		NumericVector Rllh(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat AN(RN.begin(), n, m, true);
		// n is already (nrow(data) + maxorder i.e.lag)
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
		// Unconditional intercept
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
			if(Rmodel[2]>0){
				for(i = 0; i<Rmodel[2]; i++)
				{
					AQt = AQt + Rpars[Ridx[0]+i]*( arma::trans( AZ.row(j-(i+1)) ) * AZ.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					AQt = AQt + Rpars[Ridx[2]+i]*( arma::trans( AN.row(j-(i+1)) ) * AN.row(j-(i+1)) );
				}
			}
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[1]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt/temp1;
			temp2 = arma::as_scalar( AZ.row(j) * ( arma::inv(ARt) * arma::trans(AZ.row(j))) );
			llhtemp = log(arma::det(ARt)) + temp2;
			llh+= llhtemp;
			Rllh[j] = 0.5 * llhtemp;
		}
		output[0] = QtOut;
		output[1] = Rllh;
		output[2] = 0.5 * llh;
		return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    	} catch(...) {
		::Rf_error( "rmgarch-->dccfit c++ exception (unknown reason)" );
    	}
    	return R_NilValue;
}

SEXP dccnormC2(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP Nbar, SEXP H, SEXP Z, SEXP N, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RH(H);
		Rcpp::NumericMatrix RN(N);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		// add this for safety in case something goes wrong with H dimension
		int nh = RH.nrow();
		List output(4), QtOut(n), RtOut(n);
		NumericVector Rllh(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat AH(RH.begin(), nh, m, true);
		arma::mat AN(RN.begin(), n, m, true);
		// n is already (nrow(data) + maxorder i.e.lag)
		double *Rpars = REAL(pars);
		double *Repars = REAL(epars);
		int *Rmodel = INTEGER(model);
		int *Ridx = INTEGER(idx);
		int mo = (int) Repars[2];
		double llhtemp = 0.0, llh = 0.0, temp2 = 0.0, temp3 = 0.0;
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat AUQ(m, m);
		arma::mat temp1(m, m);
		double lcons = m*log(2.0*M_PI);
		// Unconditiona/intercept
		AUQ = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;
		for(j=0; j<mo; j++)
		{
			QtOut[j] = AQbar;
			RtOut[j] = AQbar;
			Rllh[j] = 0.0;
		}
		for(j=mo; j<n; j++){
			temp1.zeros();
			temp2 = 0.0;
			temp3 = 0.0;
			llhtemp = 0.0;
			ARt.zeros();
			AQt = AUQ;
			if(Rmodel[2]>0){
				for(i = 0; i<Rmodel[2]; i++)
				{
					AQt = AQt + Rpars[Ridx[0]+i]*( arma::trans( AZ.row(j-(i+1)) ) * AZ.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					AQt = AQt + Rpars[Ridx[2]+i]*( arma::trans( AN.row(j-(i+1)) ) * AN.row(j-(i+1)) );
				}
			}
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[1]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt / temp1;
			RtOut[j] = ARt;
			temp2 = arma::as_scalar( AZ.row(j) * (arma::inv(ARt) * arma::trans(AZ.row(j)) ) );
			temp3 = 2 * log(arma::prod(AH.row(j)));
			llhtemp = lcons + temp3 + log(arma::det(ARt)) + temp2;
			llh+= llhtemp;
			Rllh[j] = -0.5 * llhtemp;
		}
		output[0] = QtOut;
		output[1] = Rllh;
		output[2] = -0.5 * llh;
		output[3] = RtOut;
		return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    	} catch(...) {
		::Rf_error( "rmgarch-->dccfit c++ exception (unknown reason)" );
    	}
    	return R_NilValue;
}


SEXP dccstudentC1(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP Nbar, SEXP Z, SEXP N, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RN(N);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		List output(3), QtOut(n);
		NumericVector Rllh(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat AN(RN.begin(), n, m, true);
		// n is already (nrow(data) + maxorder i.e.lag)
		double *Rpars = REAL(pars);
		double *Repars = REAL(epars);
		int *Rmodel = INTEGER(model);
		int *Ridx = INTEGER(idx);
		int mo = (int) Repars[2];
		double llhtemp = 0.0, llh = 0.0, temp2 = 0.0;
		double lcons = Rf_lgammafn(0.5 * (Rpars[Ridx[3]] + m)) - Rf_lgammafn(0.5 * Rpars[Ridx[3]]) - 0.5 * m * log(M_PI*(Rpars[Ridx[3]] - 2));
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat AUQ(m, m);
		arma::mat temp1(m, m);
		// Unconditional intercept
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
			if(Rmodel[2]>0){
				for(i = 0; i<Rmodel[2]; i++)
				{
					AQt = AQt + Rpars[Ridx[0]+i]*( arma::trans( AZ.row(j-(i+1)) ) * AZ.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					AQt = AQt + Rpars[Ridx[2]+i]*( arma::trans( AN.row(j-(i+1)) ) * AN.row(j-(i+1)) );
				}
			}
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[1]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt/temp1;
			temp2 = arma::as_scalar( AZ.row(j) * ( arma::inv(ARt) * arma::trans(AZ.row(j))) );
			llhtemp = lcons - 0.5 * log(arma::det(ARt)) - 0.5 * (Rpars[Ridx[3]] + m) * log(1.0 + (1.0/(Rpars[Ridx[3]] - 2)) * temp2);
			llh+= llhtemp;
			Rllh[j] = -1.0 * llhtemp;
		}
		output[0] = QtOut;
		output[1] = Rllh;
		output[2] = -1.0 * llh;
		return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    	} catch(...) {
		::Rf_error( "rmgarch-->dccfit c++ exception (unknown reason)" );
    	}
    	return R_NilValue;
}

SEXP dccstudentC2(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP Nbar, SEXP H, SEXP Z, SEXP N, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RH(H);
		Rcpp::NumericMatrix RN(N);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		// add this for safety in case something goes wrong with H dimension
		int nh = RH.nrow();
		List output(4), QtOut(n), RtOut(n);
		NumericVector Rllh(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat AH(RH.begin(), nh, m, true);
		arma::mat AN(RN.begin(), n, m, true);
		// n is already (nrow(data) + maxorder i.e.lag)
		double *Rpars = REAL(pars);
		double *Repars = REAL(epars);
		int *Rmodel = INTEGER(model);
		int *Ridx = INTEGER(idx);
		int mo = (int) Repars[2];
		double llhtemp = 0.0, llh = 0.0, temp2 = 0.0, temp3 = 0.0;
		double lcons = Rf_lgammafn(0.5 * (Rpars[Ridx[3]] + m)) - Rf_lgammafn(0.5 * Rpars[Ridx[3]]) - 0.5 * m * log(M_PI*(Rpars[Ridx[3]] - 2));
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat AUQ(m, m);
		arma::mat temp1(m, m);
		// Unconditiona/intercept
		AUQ = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;
		for(j=0; j<mo; j++)
		{
			QtOut[j] = AQbar;
			RtOut[j] = AQbar;
			Rllh[j] = 0.0;
		}
		for(j=mo; j<n; j++){
			temp1.zeros();
			temp2 = 0.0;
			temp3 = 0.0;
			llhtemp = 0.0;
			ARt.zeros();
			AQt = AUQ;
			if(Rmodel[2]>0){
				for(i = 0; i<Rmodel[2]; i++)
				{
					AQt = AQt + Rpars[Ridx[0]+i]*( arma::trans( AZ.row(j-(i+1)) ) * AZ.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					AQt = AQt + Rpars[Ridx[2]+i]*( arma::trans( AN.row(j-(i+1)) ) * AN.row(j-(i+1)) );
				}
			}
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[1]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt / temp1;
			RtOut[j] = ARt;
			temp2 = arma::as_scalar( AZ.row(j) * (arma::inv(ARt) * arma::trans(AZ.row(j))));
			temp3 = log(arma::prod(AH.row(j)));
			llhtemp = lcons - 0.5 * log(arma::det(ARt)) - temp3 - 0.5 * (Rpars[Ridx[3]] + m) * log(1.0 + (1.0/(Rpars[Ridx[3]] - 2)) * temp2);
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
		::Rf_error( "rmgarch-->dccfit c++ exception (unknown reason)" );
    	}
    	return R_NilValue;
}

SEXP dcclaplaceC1(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP Nbar, SEXP Z, SEXP N, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RN(N);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		List output(3), QtOut(n);
		NumericVector Rllh(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat AN(RN.begin(), n, m, true);
		// n is already (nrow(data) + maxorder i.e.lag)
		double *Rpars = REAL(pars);
		double *Repars = REAL(epars);
		int *Rmodel = INTEGER(model);
		int *Ridx = INTEGER(idx);
		int mo = (int) Repars[2];
		double v = 0.5 * (2 - m);
		double lcons = log(2.0) - 0.5 * m * log(2 * M_PI);
		double llhtemp = 0.0, llh = 0.0, temp2 = 0.0;
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat AUQ(m, m);
		arma::mat temp1(m, m);
		// Unconditiona/intercept
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
			if(Rmodel[2]>0){
				for(i = 0; i<Rmodel[2]; i++)
				{
					AQt = AQt + Rpars[Ridx[0]+i]*( arma::trans( AZ.row(j-(i+1)) ) * AZ.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					AQt = AQt + Rpars[Ridx[2]+i]*( arma::trans( AN.row(j-(i+1)) ) * AN.row(j-(i+1)) );
				}
			}
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[1]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt / temp1;
			temp2 = arma::as_scalar( AZ.row(j) * (arma::inv(ARt) * arma::trans(AZ.row(j))));
			llhtemp = lcons - 0.5 * log(arma::det(ARt)) + 0.5 * v * log(0.5 * temp2) + log(Rf_bessel_k(sqrt(2 * temp2), v, 1));
			llh+= llhtemp;
			Rllh[j] = -1.0 * llhtemp;
		}
		output[0] = QtOut;
		output[1] = Rllh;
		output[2] = -1.0 * llh;
		return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    	} catch(...) {
		::Rf_error( "rmgarch-->dccfit c++ exception (unknown reason)" );
    	}
    	return R_NilValue;
}


SEXP dcclaplaceC2(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP Nbar, SEXP H, SEXP Z, SEXP N, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RH(H);
		Rcpp::NumericMatrix RN(N);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		// add this for safety in case something goes wrong with H dimension
		int nh = RH.nrow();
		List output(4), QtOut(n), RtOut(n);
		NumericVector Rllh(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat AH(RH.begin(), nh, m, true);
		arma::mat AN(RN.begin(), n, m, true);
		// n is already (nrow(data) + maxorder i.e.lag)
		double *Rpars = REAL(pars);
		double *Repars = REAL(epars);
		int *Rmodel = INTEGER(model);
		int *Ridx = INTEGER(idx);
		int mo = (int) Repars[2];
		double llhtemp = 0.0, llh = 0.0, temp2 = 0.0, temp3 = 0.0;
		double v = 0.5 * (2 - m);
		double lcons = log(2.0) - 0.5 * m * log(2 * M_PI);
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat AUQ(m, m);
		arma::mat temp1(m, m);
		// Unconditiona/intercept
		AUQ = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;
		// rem: newmat matrices start at one not zero
		// populate Dt/Dtlist
		for(j=0; j<mo; j++)
		{
			QtOut[j] = AQbar;
			RtOut[j] = AQbar;
			Rllh[j] = 0.0;
		}
		for(j=mo; j<n; j++){
			temp1.zeros();
			temp2 = 0.0;
			temp3 = 0.0;
			llhtemp = 0.0;
			ARt.zeros();
			AQt = AUQ;
			if(Rmodel[2]>0){
				for(i = 0; i<Rmodel[2]; i++)
				{
					AQt = AQt + Rpars[Ridx[0]+i]*( arma::trans( AZ.row(j-(i+1)) ) * AZ.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					AQt = AQt + Rpars[Ridx[2]+i]*( arma::trans( AN.row(j-(i+1)) ) * AN.row(j-(i+1)) );
				}
			}
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[1]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt / temp1;
			RtOut[j] = ARt;
			temp2 = arma::as_scalar( AZ.row(j) * (arma::inv(ARt) * arma::trans(AZ.row(j))));
			temp3 = log(arma::prod(AH.row(j)));
			llhtemp = lcons - 0.5 * log(arma::det(ARt)) - temp3 + 0.5 * v * log(0.5 * temp2) + log(Rf_bessel_k(sqrt(2 * temp2), v, 1));
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
		::Rf_error( "rmgarch-->dccfit c++ exception (unknown reason)" );
    	}
    	return R_NilValue;
}

SEXP dccfilter(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP Rbar, SEXP Nbar, SEXP Z, SEXP epars)
{
	try {
		Rcpp::NumericMatrix RQbar(Qbar);
		Rcpp::NumericMatrix RNbar(Nbar);
		Rcpp::NumericMatrix RRbar(Rbar);
		Rcpp::NumericMatrix RZ(Z);
		int m = RQbar.nrow(), n = RZ.nrow(), i, j;
		List output(3), RtOut(n), QtOut(n), ZtOut(n);
		NumericVector Rllh(n);
		arma::mat AQbar(RQbar.begin(), m, m, true);
		arma::mat ANbar(RNbar.begin(), m, m, true);
		arma::mat ARbar(RRbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
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
		for(j=0; j<mo; j++)
		{
			QtOut[j] = AQbar;
			RtOut[j] = ARbar;
			AZt.row(j) = AZ.row(j);
			ANt.row(j) = ( 0.5*(AZt.row(j)/abs(AZt.row(j))-1) )%AZt.row(j);
		}
		for(j=mo; j<n; j++)
		{
			AQt = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;
			if(Rmodel[2]>0){
				for(i = 0; i<Rmodel[2]; i++)
				{
					AQt = AQt + Rpars[Ridx[0]+i]*( arma::trans( AZt.row(j-(i+1)) ) * AZt.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					AQt = AQt + Rpars[Ridx[2]+i]*( arma::trans( ANt.row(j-(i+1)) ) * ANt.row(j-(i+1)) );
				}
			}
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[1]+i] * Qtemp;
				}
			}
			QtOut[j] =AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt / temp1;
			RtOut[j] = ARt;
			eig_sym(eigval, eigvec, ARt);
			temp1 = ( eigvec * arma::diagmat( arma::sqrt( eigval ) ) * arma::inv( eigvec ) );
			AZt.row(j) = AZ.row(j) * temp1;
			// To obtain the negative of the Zs sign(z) = (z/abs(z)-1)/2
			ANt.row(j) = ( 0.5*(AZt.row(j)/abs(AZt.row(j))-1) )%AZt.row(j);
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

SEXP dccsimmvn(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP preQ, SEXP Rbar, SEXP Nbar, SEXP Z, SEXP NZ, SEXP epars)
{
	// Z are the pre-Z's (distribution specific)
	// NZ are the normal Z's used by all distributions
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
		arma::mat ApreQ(RpreQ.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat ANZ(RNZ.begin(), n, m, true);
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
		arma::rowvec tmpz(m);
		for(j=0; j<mo; j++)
		{
			QtOut[j] = ApreQ;
			RtOut[j] = ARbar;
			AZt.row(j) = AZ.row(j);
			ANt.row(j) = ( 0.5*(AZt.row(j)/abs(AZt.row(j)+1E-12)-1) )%AZt.row(j);
		}
		for(j=mo; j<n; j++)
		{
			AQt = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;
			if(Rmodel[2]>0){
				for(i = 0; i<Rmodel[2]; i++)
				{
					AQt = AQt + Rpars[Ridx[0]+i]*( arma::trans( AZt.row(j-(i+1)) ) * AZt.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					AQt = AQt + Rpars[Ridx[2]+i]*( arma::trans( ANt.row(j-(i+1)) ) * ANt.row(j-(i+1)) );
				}
			}
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[1]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt / temp1;
			arma::mat Rtmp(m,m);
			// enforce this in case of problems
			Rtmp = arma::symmatu(ARt);
			RtOut[j] = Rtmp;
			//eig_sym(eigval, eigvec, Rtmp);
			// generate conditional standard random deviates using the multivariate normal method
			// This will only be an approximation for all other distributions used.
			//temp1 = ( eigvec * arma::diagmat( arma::sqrt( eigval ) ) * arma::inv( eigvec ) );
			//AZt.row(j) = ARes.row(j) * temp1;
			tmpz = ANZ.row(j);
			AZt.row(j) = rmvnormx(Rtmp, tmpz);
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

SEXP dccsimmvl(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP preQ, SEXP Rbar, SEXP Nbar, SEXP Z, SEXP NZ, SEXP epars)
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
		arma::mat ApreQ(RpreQ.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat ANZ(RNZ.begin(), n, m, true);
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
		arma::rowvec tmpz(m);
		for(j=0; j<mo; j++)
		{
			QtOut[j] = ApreQ;
			RtOut[j] = ARbar;
			AZt.row(j) = AZ.row(j);
			ANt.row(j) = ( 0.5*(AZt.row(j)/abs(AZt.row(j)+1E-12)-1) )%AZt.row(j);
		}
		for(j=mo; j<n; j++)
		{
			AQt = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;
			if(Rmodel[2]>0){
				for(i = 0; i<Rmodel[2]; i++)
				{
					AQt = AQt + Rpars[Ridx[0]+i]*( arma::trans( AZt.row(j-(i+1)) ) * AZt.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					AQt = AQt + Rpars[Ridx[2]+i]*( arma::trans( ANt.row(j-(i+1)) ) * ANt.row(j-(i+1)) );
				}
			}
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[1]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt / temp1;
			arma::mat Rtmp(m,m);
			// enforce this in case of problems
			Rtmp = arma::symmatu(ARt);
			RtOut[j] = Rtmp;
			//eig_sym(eigval, eigvec, Rtmp);
			// generate conditional standard random deviates using the multivariate normal method
			// This will only be an approximation for all other distributions used.
			//temp1 = ( eigvec * arma::diagmat( arma::sqrt( eigval ) ) * arma::inv( eigvec ) );
			//AZt.row(j) = ARes.row(j) * temp1;
			tmpz = ANZ.row(j);
			AZt.row(j) = rmvlx(Rtmp, tmpz);
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

SEXP dccsimmvt(SEXP model, SEXP pars, SEXP idx, SEXP Qbar, SEXP preQ, SEXP Rbar, SEXP Nbar, SEXP Z, SEXP NZ, SEXP epars)
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
		arma::mat ApreQ(RpreQ.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat ANZ(RNZ.begin(), n, m, true);
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
		double nu = Rpars[Ridx[3]];
		arma::rowvec tmpz(m);
		for(j=0; j<mo; j++)
		{
			QtOut[j] = ApreQ;
			RtOut[j] = ARbar;
			AZt.row(j) = AZ.row(j);
			ANt.row(j) = ( 0.5*(AZt.row(j)/abs(AZt.row(j)+1E-12)-1) )%AZt.row(j);
		}
		for(j=mo; j<n; j++)
		{
			AQt = AQbar * (1.0 - Repars[0]) - Repars[1] * ANbar;
			if(Rmodel[2]>0){
				for(i = 0; i<Rmodel[2]; i++)
				{
					AQt = AQt + Rpars[Ridx[0]+i]*( arma::trans( AZt.row(j-(i+1)) ) * AZt.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					AQt = AQt + Rpars[Ridx[2]+i]*( arma::trans( ANt.row(j-(i+1)) ) * ANt.row(j-(i+1)) );
				}
			}
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + Rpars[Ridx[1]+i] * Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt / temp1;
			arma::mat Rtmp(m,m);
			// enforce this in case of problems
			Rtmp = arma::symmatu(ARt);
			RtOut[j] = Rtmp;
			//eig_sym(eigval, eigvec, Rtmp);
			// generate conditional standard random deviates using the multivariate normal method
			// This will only be an approximation for all other distributions used.
			//temp1 = ( eigvec * arma::diagmat( arma::sqrt( eigval ) ) * arma::inv( eigvec ) );
			//AZt.row(j) = ARes.row(j) * temp1;
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

SEXP fdccnormC1(SEXP model, SEXP A, SEXP B, SEXP C, SEXP Qbar, SEXP Z)
{
	try {
			Rcpp::NumericMatrix RA(A);
			Rcpp::NumericMatrix RB(B);
			Rcpp::NumericMatrix RC(C);
			Rcpp::NumericMatrix RZ(Z);
			Rcpp::NumericMatrix RQ(Qbar);
			int m = RZ.ncol(), n = RZ.nrow(), i, j;
			List output(3), QtOut(n);
			NumericVector Rllh(n);
			int rA = RA.ncol(), rB = RB.ncol(), rC = RC.nrow();
			arma::mat AA(RA.begin(), m, rA, true);
			arma::mat AB(RB.begin(), m, rB, true);
			arma::mat AC(RC.begin(), m, rC, true);
			arma::mat AZ(RZ.begin(), n, m, true);
			arma::mat AQ(RQ.begin(), m, m, true);
			// n is already (nrow(data) + maxorder i.e.lag)
			int *Rmodel = INTEGER(model);
			int s = Rmodel[3] > Rmodel[4] ? Rmodel[3] : Rmodel[4];
			double llhtemp = 0.0, llh = 0.0, temp2 = 0.0;
			// empty matrices
			arma::mat AQt(m, m);
			arma::mat ARt(m, m);
			arma::mat temp1(m, m);
			arma::mat Qtemp(m,m);
			// Intercept
			arma::mat Q = AC%AQ;
			for(j=0; j<s; j++)
			{
				QtOut[j] = AQ;
				Rllh[j] = 0.0;
			}
			for(j=s; j<n; j++)
			{
				temp2 = 0.0;
				llhtemp = 0.0;
				temp1.zeros();
				ARt.zeros();
				AQt = Q;
				if(Rmodel[3]>0){
					for(i = 0; i<Rmodel[3]; i++)
					{
						AQt = AQt + (AA.col(i)*arma::trans(AA.col(i))) % ( arma::trans( AZ.row(j-(i+1)) ) * AZ.row(j-(i+1)) );
					}
				}
				if(Rmodel[4]>0){
					for(i = 0; i<Rmodel[4]; i++)
					{
						Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
						AQt = AQt + (AB.col(i)*arma::trans(AB.col(i))) % Qtemp;
					}
				}
				QtOut[j] = AQt;
				temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
				ARt = AQt/temp1;
				temp2 = arma::as_scalar( AZ.row(j) * ( arma::inv(ARt) * arma::trans(AZ.row(j))) );
				llhtemp = log(arma::det(ARt)) + temp2;
				llh+= llhtemp;
				Rllh[j] = 0.5 * llhtemp;
			}
			output[0] = QtOut;
			output[1] = Rllh;
			output[2] = 0.5 * llh;
			return(output);
		} catch( std::exception &ex ) {
			forward_exception_to_r( ex );
	    	} catch(...) {
			::Rf_error( "rmgarch-->fdccfit c++ exception (unknown reason)" );
	    	}
	    	return R_NilValue;
}

SEXP fdccnormC2(SEXP model, SEXP A, SEXP B, SEXP C, SEXP Qbar, SEXP H, SEXP Z)
{
	try {
		Rcpp::NumericMatrix RA(A);
		Rcpp::NumericMatrix RB(B);
		Rcpp::NumericMatrix RC(C);
		Rcpp::NumericMatrix RH(H);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RQ(Qbar);
		int m = RZ.ncol(), n = RZ.nrow(), i, j;
		int nh = RH.nrow();
		List output(4), QtOut(n), RtOut(n);
		int rA = RA.ncol(), rB = RB.ncol(), rC = RC.ncol();
		NumericVector Rllh(n);
		arma::mat AA(RA.begin(), m, rA, true);
		arma::mat AB(RB.begin(), m, rB, true);
		arma::mat AC(RC.begin(), m, rC, true);
		arma::mat AQ(RQ.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat AH(RH.begin(), nh, m, true);
		// n is already (nrow(data) + maxorder i.e.lag)
		int *Rmodel = INTEGER(model);
		int s = Rmodel[3] > Rmodel[4] ? Rmodel[3] : Rmodel[4];
		double llhtemp = 0.0, llh = 0.0, temp2 = 0.0, temp3 = 0.0;
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat temp1(m, m);
		arma::mat Qtemp(m,m);
		// Intercept
		arma::mat Q = AC%AQ;
		double lcons = m*log(2.0*M_PI);
		for(j=0; j<s; j++)
		{
			QtOut[j] = AQ;
			RtOut[j] = AQ;
			Rllh[j] = 0.0;
		}
		for(j=s; j<n; j++){
			temp1.zeros();
			temp2 = 0.0;
			temp3 = 0.0;
			llhtemp = 0.0;
			ARt.zeros();
			AQt = Q;
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					AQt = AQt +  (AA.col(i)*arma::trans(AA.col(i))) % ( arma::trans( AZ.row(j-(i+1)) ) * AZ.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt +  (AB.col(i)*arma::trans(AB.col(i))) % Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt / temp1;
			RtOut[j] = ARt;
			temp2 = arma::as_scalar( AZ.row(j) * (arma::inv(ARt) * arma::trans(AZ.row(j)) ) );
			temp3 = 2 * log(arma::prod(AH.row(j)));
			llhtemp = lcons + temp3 + log(arma::det(ARt)) + temp2;
			llh+= llhtemp;
			Rllh[j] = -0.5 * llhtemp;
		}
		output[0] = QtOut;
		output[1] = Rllh;
		output[2] = -0.5 * llh;
		output[3] = RtOut;
		return(output);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
    	} catch(...) {
		::Rf_error( "rmgarch-->fdccfit c++ exception (unknown reason)" );
    	}
    	return R_NilValue;
}

SEXP fdccsimmvn(SEXP model, SEXP A, SEXP B, SEXP C, SEXP Qbar, SEXP preQ, SEXP Rbar, SEXP Z, SEXP NZ)
{
	// Z are the pre-Z's (distribution specific)
	// NZ are the normal Z's used by all distributions
	try {
		Rcpp::NumericMatrix RA(A);
		Rcpp::NumericMatrix RB(B);
		Rcpp::NumericMatrix RC(C);
		Rcpp::NumericMatrix RRbar(Rbar);
		Rcpp::NumericMatrix RNZ(NZ);
		Rcpp::NumericMatrix RZ(Z);
		Rcpp::NumericMatrix RQ(Qbar);
		Rcpp::NumericMatrix RPREQ(preQ);
		int m = RRbar.nrow(), n = RZ.nrow(), i, j;
		int rA = RA.ncol(), rB = RB.ncol(), rC = RC.ncol();
		List output(3), RtOut(n), QtOut(n), ZtOut(n);
		arma::mat AA(RA.begin(), m, rA, true);
		arma::mat AB(RB.begin(), m, rB, true);
		arma::mat AC(RC.begin(), m, rC, true);
		arma::mat ARbar(RRbar.begin(), m, m, true);
		arma::mat AZ(RZ.begin(), n, m, true);
		arma::mat ANZ(RNZ.begin(), n, m, true);
		arma::mat AQ(RQ.begin(), m, m, true);
		arma::mat APREQ(RPREQ.begin(), m, m, true);
		int *Rmodel = INTEGER(model);
		int s = Rmodel[3] > Rmodel[4] ? Rmodel[3] : Rmodel[4];
		// empty matrices
		arma::mat AQt(m, m);
		arma::mat ARt(m, m);
		arma::mat temp1(m, m);
		arma::mat AZt(n, m);
		arma::vec eigval(m);
		arma::mat eigvec(m, m);
		arma::rowvec tmpz(m);
		arma::mat Q = AC%AQ;
		for(j=0; j<s; j++)
		{
			QtOut[j] = APREQ;
			RtOut[j] = ARbar;
			AZt.row(j) = AZ.row(j);
		}
		for(j=s; j<n; j++)
		{
			AQt = Q;
			if(Rmodel[3]>0){
				for(i = 0; i<Rmodel[3]; i++)
				{
					AQt = AQt + (AA.col(i)*arma::trans(AA.col(i))) % ( arma::trans( AZt.row(j-(i+1)) ) * AZt.row(j-(i+1)) );
				}
			}
			if(Rmodel[4]>0){
				for(i = 0; i<Rmodel[4]; i++)
				{
					arma::mat Qtemp = Rcpp::as<arma::mat>(QtOut[j-(i+1)]);
					AQt = AQt + (AB.col(i)*arma::trans(AB.col(i))) % Qtemp;
				}
			}
			QtOut[j] = AQt;
			temp1 = arma::sqrt(AQt.diag()) * arma::trans(arma::sqrt(AQt.diag()));
			ARt = AQt / temp1;
			arma::mat Rtmp(m,m);
			// enforce this in case of problems
			Rtmp = arma::symmatu(ARt);
			RtOut[j] = Rtmp;
			tmpz = ANZ.row(j);
			AZt.row(j) = rmvnormx(Rtmp, tmpz);
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
