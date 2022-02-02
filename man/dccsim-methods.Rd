\name{dccsim-methods}
\docType{methods}
\alias{dccsim}
\alias{dccsim,ANY-method}
\alias{dccsim,DCCspec-method}
\alias{dccsim,DCCfit-method}
\title{function: DCC-GARCH Simulation}
\description{
Method for creating a DCC-GARCH simulation object.
}
\usage{
dccsim(fitORspec, n.sim = 1000, n.start = 0, m.sim = 1, 
startMethod = c("unconditional", "sample"), presigma = NULL, preresiduals = NULL, 
prereturns = NULL, preQ = NULL, preZ = NULL, Qbar = NULL, Nbar = NULL, 
rseed = NULL, mexsimdata = NULL, vexsimdata = NULL, cluster = NULL, 
VAR.fit = NULL, prerealized = NULL, ...)
}
\arguments{
  	\item{fitORspec}{ A \code{\linkS4class{DCCspec}} or \code{\linkS4class{DCCfit}} 
  	object created by calling either \code{\link{dccspec}} with fixed parameters 
  	or \code{\link{dccfit}}.}
	\item{n.sim}{ The simulation horizon.}
	\item{n.start}{ The burn-in sample.}
	\item{m.sim}{ The number of simulations.}
	\item{startMethod}{ Starting values for the simulation. Valid methods are 
	\dQuote{unconditional} for the expected values given the density, and 
	\dQuote{sample} for the ending values of the actual data from the fit 
	object (for the dispatch method using a specification, \dQuote{sample} is 
	not relevant).}
	\item{presigma}{ Allows the starting sigma values to be provided by the user
	for the univariate GARCH dynamics.}
	\item{prereturns}{ Allows the starting return data to be provided by the 
	user for the conditional mean simulation.}
	\item{preresiduals}{ Allows the starting residuals to be provided by the 
	user and used in the GARCH dynamics simulation.}
	\item{preQ}{ Allows the starting \sQuote{DCC-Q} value to be provided by the 
	user and though unnecessary for the first 1-ahead simulation using the 
	\dQuote{sample} option in the \code{startMethod}, this is key to obtaining
	a rolling n-ahead forecast type simulation (see details below).}
	\item{preZ}{ Allows the starting standardized residuals to be provided by 
	the user and though unnecessary for the first 1-ahead simulation using 
	the \dQuote{sample} option in the \code{startMethod}, this is key to 
	obtaining a rolling n-ahead forecast type simulation (see details below).}
	\item{Qbar}{ The DCC dynamics unconditional Q matrix, required for the 
	specification dispatch method.}
	\item{Nbar}{ The aDCC dynamics unconditional asymmetry matrix, required for 
	the specification dispatch method.}
	\item{rseed}{ Optional seeding value(s) for the random number generator. 
	For m.sim>1, it is possible to provide either a single seed to initialize 
	all values, or one seed per separate simulation (i.e. m.sim seeds). However, 
	in the latter case this may result in some slight overhead depending on
	how large m.sim is.}
	\item{mexsimdata}{ A list (equal to the number of asset) of matrices of 
	simulated external regressor-in-mean data with row length equal to 
	n.sim + n.start. If the fit object contains external regressors in the mean 
	equation, this must be provided else will be assumed to be zero.}
	\item{vexsimdata}{ A list (equal to the number of asset) of matrices of 
	simulated external regressor-in-variance data with row length equal to 
	n.sim + n.start. If the fit object contains external regressors in the 
	variance equation, this must be provided else will be assumed to be zero.}
  	\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
  	the parallel package. If it is not NULL, then this will be used for parallel 
  	estimation (remember to stop the cluster on completion).}
	\item{VAR.fit}{ An VAR.fit list returned from calling the 
	\code{\link{varxfilter}} or \code{\link{varxfit}} function with \code{postpad} 
	set to \dQuote{constant}. This is required for the specification dispatch
	method.}
	\item{prerealized}{Allows the starting realized volatility values to be provided 
	by the user for the univariate GARCH dynamics.}
  	\item{...}{.}
}
\value{
A \code{\linkS4class{DCCsim}} object containing details of the DCC-GARCH 
simulation.
}
\details{
In order to pass a correct specification to the filter routine, you must ensure 
that it contains the appropriate \sQuote{fixed.pars} in both the multivariate 
DCC part of the specification as well as the multiple univariate specification 
part, for which the method \code{\link{setfixed<-}} should be used.
}
\author{Alexios Galanos}
\keyword{methods}
