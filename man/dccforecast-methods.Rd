\name{dccforecast-methods}
\docType{methods}
\alias{dccforecast}
\alias{dccforecast,ANY-method}
\alias{dccforecast,DCCfit-method}
\title{function: DCC-GARCH Forecast}
\description{
Method for creating a DCC-GARCH forecast object.
}
\usage{
dccforecast(fit, n.ahead = 1, n.roll = 0, 
external.forecasts = list(mregfor = NULL, vregfor = NULL), cluster = NULL, ...) 
}
\arguments{
  	\item{fit}{ A \code{\linkS4class{DCCfit}} object created by calling 
  	\code{\link{dccfit}}.}
  	\item{n.ahead}{ The forecast horizon.}
  	\item{n.roll}{ The no. of rolling forecasts to create beyond the first one 
  	(see details).}
  	\item{external.forecasts}{ A list with forecasts for the external regressors 
  	in the mean and/or variance equations if specified (see details).}
  	\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
  	the parallel package. If it is not NULL, then this will be used for parallel 
  	estimation (remember to stop the cluster on completion).}
  	\item{...}{.}
}
\value{
A \code{\linkS4class{DCCforecast}} object containing details of the DCC-GARCH 
forecast.
}
\details{
When using \code{n.roll}, it is assumed that \code{\link{dccfit}} was called 
with argument \sQuote{out.sample} being large enough to cover n-rolling 
forecasts.\cr 
When n.roll = 0, all forecasts are based on an unconditional n-ahead forecast 
routine based on the approximation method described in ENGLE and SHEPPARD (2001) 
paper (see reference below). If any external regressors are present, then the 
user must pass in their unconditional forecasts in the \sQuote{external.forecasts} 
list, as matrices with dimensions equal to n.ahead x n.assets. This assumes
that the univariate GARCH specifications share common external regressors 
(this may change in the future).\cr
When n.roll>0 and n.ahead = 1, then this is a pure rolling forecast based on the 
available out.sample data provided for in the call to the fit routine. It is 
also assumed that if any external regressors were passed to the fit routine that 
they contained enough values to cover the out.sample period so that they could 
be used in this forecast scenario.\cr 
The case of n.roll > 0 AND n.ahead > 1 is not implemented.\cr
}
\references{
Engle, R.F. and Sheppard, K. 2001, Theoretical and empirical properties of 
dynamic conditional correlation multivariate GARCH, \emph{NBER Working Paper}.\cr
}
\author{Alexios Galanos}
\keyword{methods}
