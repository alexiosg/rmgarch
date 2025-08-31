\name{dccfit-methods}
\docType{methods}
\alias{dccfit}
\alias{dccfit,ANY-method}
\alias{dccfit,DCCspec-method}
\title{function: DCC-GARCH Fit}
\description{
Method for creating a DCC-GARCH fit object.
}
\usage{
dccfit(spec, data, out.sample = 0, solver = "solnp", solver.control = list(),
fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE),
cluster = NULL, fit = NULL, VAR.fit = NULL, realizedVol = NULL, ...)
}
\arguments{
  	\item{spec}{ A \code{\linkS4class{DCCspec}} object created by calling
  	\code{\link{dccspec}}.}
  	\item{data}{ A multivariate data object of class xts or one which can be
  	coerced to such.}
  	\item{out.sample}{ A positive integer indicating the number of periods
  	before the last to keep for out of sample forecasting.}
  	\item{solver}{ Either \dQuote{nlminb}, \dQuote{solnp}, \dQuote{gosolnp} or
  	\dQuote{lbfgs}. It can also optionally be a vector of length 2 with the first
  	solver being used for the first stage univariate GARCH estimation (in which
  	case the option of \dQuote{hybrid} is also available).}
 	\item{solver.control}{Control arguments list passed to optimizer.}
	\item{fit.control}{Control arguments passed to the fitting routine.
	The \sQuote{eval.se} option determines whether standard errors are calculated
	(see details below). The \sQuote{stationarity} option is for the univariate
	stage GARCH fitting routine, whilst for the second stage DCC this is
	imposed by design. The \sQuote{scale} option is also for the first stage
	univariate GARCH fitting routine.}
  	\item{cluster}{ A cluster object created by calling \code{makeCluster} from
  	the parallel package. If it is not NULL, then this will be used for parallel
  	estimation (remember to stop the cluster on completion).}
  \item{fit}{ (optional) A previously estimated univariate
  \code{uGARCHmultifit} object from the \pkg{rugarch} package (see details).}
	\item{VAR.fit}{ (optional) A previously estimated VAR object returned from
	calling the \code{\link{varxfit}} function.}
	\item{realizedVol}{ Required xts matrix for the realGARCH model.}
	\item{...}{ . }
}
\value{
A \code{\linkS4class{DCCfit}} object containing details of the DCC-GARCH fit.
}
\author{Alexios Galanos}
\details{
The 2-step DCC estimation fits a GARCH-Normal model to the univariate data and
then proceeds to estimate the second step based on the chosen multivariate
distribution. Because of this 2-step approach, standard errors are expensive to
calculate and therefore the use of parallel functionality, built into both the
fitting and standard error calculation routines is key. The switch to turn off
the calculation of standard errors through the \sQuote{fit.control} option could
be quite useful in rolling estimation such as in the \code{\link{dccroll}}
routine.\cr
The optional \sQuote{fit} argument allows to pass your own
\code{uGARCHmultifit} object from the \pkg{rugarch} package instead of having the
routine estimate it. This is very useful in cases
of multiple use of the same fit and problems in convergence which might require
a more hands on approach to the univariate fitting stage. However, it is up to
the user to ensure consistency between the \sQuote{fit} and supplied \sQuote{spec}.
}
\note{
There is no check on the VAR.fit list passed to the method so particular care
should be exercised so that the same data used in the fitting routine is also
used in the VAR fit routine. This this must have been called with the option
\code{postpad} \sQuote{constant}. The ability to pass this list of the
pre-calculated VAR model is particularly useful when comparing different models
(such as copula-GARCH, GO-GARCH etc) using the same dataset and VAR method (i.e.
the same first stage conditional mean filtration). Though the classical VAR
estimation is very fast and may not require this extra step, the robust method
is slow and therefore benefits from calculating this only once.\cr
For extensive examples look in the \sQuote{rmgarch.tests} folder.
}
\keyword{methods}
