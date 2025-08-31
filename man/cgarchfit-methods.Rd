\name{cgarchfit-methods}
\docType{methods}
\alias{cgarchfit}
\alias{cgarchfit,ANY-method}
\alias{cgarchfit,cGARCHspec-method}
\title{function: Copula-GARCH Fit}
\description{
Method for creating a Copula-GARCH fit object.
}
\usage{
cgarchfit(spec, data, spd.control = list(lower = 0.1, upper = 0.9, type = "pwm",
kernel = "epanech"), fit.control = list(eval.se = TRUE, stationarity = TRUE,
scale = FALSE), solver = "solnp", solver.control = list(), out.sample = 0,
cluster = NULL, fit = NULL, VAR.fit = NULL, realizedVol = NULL,...)
}
\arguments{
  	\item{spec}{ A \code{\linkS4class{cGARCHspec}} A \code{\linkS4class{cGARCHspec}}
  	object created by calling \code{\link{cgarchspec}}.}
  	\item{data}{ A multivariate xts data object or one which can be coerced to
  	such.}
  	\item{out.sample}{ A positive integer indicating the number of periods before
  	the last to keep for out of sample forecasting.}
  	\item{solver}{ Either \dQuote{nlminb}, \dQuote{solnp}, \dQuote{gosolnp} or
  	\dQuote{lbfgs}. It can also optionally be a vector of length 2 with the first
  	solver being used for the first stage univariate GARCH estimation (in which
  	case the option of \dQuote{hybrid} is also available).}
 	\item{solver.control}{Control arguments list passed to optimizer.}
	\item{fit.control}{Control arguments passed to the fitting routine.
	The \sQuote{eval.se} option determines whether standard errors are calculated
	(see details below). The \sQuote{scale} option is for the first stage
	univariate GARCH fitting routine.}
  	\item{cluster}{ A cluster object created by calling \code{makeCluster} from
  	the parallel package. If it is not NULL, then this will be used for parallel
  	estimation (remember to stop the cluster on completion).}
  \item{fit}{ (optional) A previously estimated univariate
  \code{uGARCHmultifit} object from the \pkg{rugarch} package (see details).}
	\item{VAR.fit}{ (optional) A previously estimated VAR list returned from
	calling the \code{\link{varxfit}} function.}
	\item{spd.control}{ If the spd transformation was chosen in the
	specification, the spd.control passes its arguments to the
	\code{spdfit} routine of the \code{spd} package.}
	\item{realizedVol}{ Required xts matrix for the realGARCH model.}
	\item{...}{ . }
}
\value{
A \code{\linkS4class{cGARCHfit}} Object containing details of the Copula-GARCH
fit.
}
\author{Alexios Galanos}
\details{
The Copula-GARCH models implemented can either be time-varying of DCC variety
else static. The multivariate Normal and Student distributions are used in the
construction of the copulas, and 3 transformation methods are available
(parametric, semi-parametric, and empirical). For the semi-parametric case the
\sQuote{spd} package of the author is available to download from CRAN and fits a
Gaussian kernel in the interior and gpd distribution for the tails (see that
package for more details).\cr
The static copula allows for the estimation of the correlation matrix either by
Maximum Likelihood or the Kendall method for the multivariate Student.\cr
Note that the \sQuote{cgarchfit} method will assign to the global environment
the \code{uGARCHmultifit} object from the \pkg{rugarch} package once that is
estimated in order to allow the routine to be restarted should something go
wrong (it should show up as \sQuote{.fitlist}).
}
\note{
There is no check on the VAR.fit list passed to the method so particular care
should be exercised so that the same data used in the fitting routine is also
used in the VAR fit routine. This must have been called with the option
\code{postpad} \sQuote{constant}. The ability to pass this list of the
pre-calculated VAR model is particularly useful when comparing different models
(such as DCC GARCH, GO GARCH etc) using the same dataset and VAR method (i.e.
the same first stage conditional mean filtration). Though the classical VAR
estimation is very fast and may not require this extra step, the robust method
is slow and therefore benefits from calculating this only once.\cr
For extensive examples look in the \sQuote{rmgarch.tests} folder.
}
\keyword{methods}
