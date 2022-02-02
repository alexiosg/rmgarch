\name{cgarchsim-methods}
\docType{methods}
\alias{cgarchsim}
\alias{cgarchsim,ANY-method}
\alias{cgarchsim,cGARCHfit-method}
\title{function: Copula-GARCH Simulation}
\description{
Method for creating a Copula-GARCH simulation object.
}
\usage{
cgarchsim(fit, n.sim = 1000, n.start = 0, m.sim = 1, 
startMethod = c("unconditional", "sample"), presigma = NULL, preresiduals = NULL, 
prereturns = NULL, preR = NULL, preQ = NULL, preZ = NULL, rseed = NULL, 
mexsimdata = NULL, vexsimdata = NULL, cluster = NULL, only.density = FALSE, 
prerealized = NULL, ...)
}
\arguments{
  	\item{fit}{ 
  	A \code{\linkS4class{cGARCHfit}} object created by calling \code{\link{cgarchfit}}.}
	\item{n.sim}{ The simulation horizon.}
	\item{n.start}{ The burn-in sample.}
	\item{m.sim}{ The number of simulations.}
	\item{startMethod}{ Starting values for the simulation. Valid methods are 
	\sQuote{unconditional} for the expected values given the density, 
	and \sQuote{sample} for the ending values of the actual data from the fit 
	object. This is mostly related to the univariate GARCH dynamics.}
	\item{presigma}{ Allows the starting sigma values to be provided by the user
	for the univariate GARCH dynamics.}
	\item{prereturns}{ Allows the starting return data to be provided by the 
	user for the conditional mean simulation.}
	\item{preresiduals}{ Allows the starting residuals to be provided by the 
	user and used in the GARCH dynamics simulation.}
	\item{preR}{ Allows the starting correlation to be provided by the user and
	mostly useful for the static copula.}
	\item{preQ}{ Allows the starting \sQuote{DCC-Q} value to be provided by the 
	user and though unnecessary for the first 1-ahead simulation using the 
	\dQuote{sample} option in the \code{startMethod}, this is key to obtaining
	a rolling n-ahead forecast type simulation (see details below).}
	\item{preZ}{ Allows the starting transformed standardized residuals (used in
	the DCC model) to be provided by the user and though unnecessary for the 
	first 1-ahead simulation using the \dQuote{sample} option in the 
	\code{startMethod}, this is key to obtaining a rolling n-ahead forecast type 
	simulation (see details below).}
	\item{rseed}{ Optional seeding value(s) for the random number generator. 
	This should be of length equal to m.sim.}
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
  	\item{only.density}{Whether to return only the simulated returns (discrete
  	time approximation to the multivariate density). This is sometimes useful
  	in order to control memory management for large simulations not requiring
  	any other information.}
  	\item{prerealized}{Allows the starting realized volatility values to be provided 
	by the user for the univariate GARCH dynamics.}
  	\item{...}{.}
}
\value{
A \code{\linkS4class{cGARCHsim}} object containing details of the Copula-GARCH 
simulation.
}
\details{
Since there is no explicit forecasting routine, the user should use this method 
for incrementally building up n-ahead forecasts by simulating 1-ahead, obtaining 
the means of the returns, sigma, Rho etc and feeding them to the next round of 
simulation as starting values. The  \sQuote{rmgarch.tests} folder contains 
specific examples which illustrate this particular point.
}
\references{
Joe, H. \emph{Multivariate Models and Dependence Concepts}, 1997, 
Chapman \& Hall, London.\cr
Genest, C., Ghoudi, K. and Rivest, L. \emph{A semiparametric estimation 
procedure of dependence parameters in multivariate families of distributions}, 
1995, Biometrika, 82, 543-552.\cr
}
\author{Alexios Galanos}
\keyword{methods}
