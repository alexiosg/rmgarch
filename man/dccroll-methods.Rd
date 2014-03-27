\name{dccroll-methods}
\docType{methods}
\alias{dccroll}
\alias{dccroll,ANY-method}
\alias{dccroll,DCCspec-method}
\title{function: DCC-GARCH Rolling Forecast}
\description{
Method for creating a DCC-GARCH rolling forecast object.
}
\usage{
dccroll(spec, data, n.ahead = 1, forecast.length = 50, refit.every = 25, 
n.start = NULL, refit.window = c("recursive", "moving"), window.size = NULL, 
solver = "solnp", solver.control = list(), 
fit.control = list(eval.se = TRUE, stationarity = TRUE, scale = FALSE), 
cluster = NULL, save.fit = FALSE, save.wdir = NULL, realizedVol = NULL, ...)
}
\arguments{
	\item{spec}{ A \code{\linkS4class{DCCspec}} object with fixed parameters.}
	\item{data}{ A multivariate xts dataset or one which can be coerced to such.}
	\item{n.ahead}{ The number of periods to forecast. }
	\item{forecast.length}{ The length of the total forecast for which out of 
	sample data from the dataset will be used for testing.}
	\item{n.start}{ Instead of forecast.length, this determines the starting 
	point in the dataset from which to initialize the rolling forecast.}
  	\item{refit.every}{ Determines every how many periods the model is 
  	re-estimated.}
	\item{refit.window}{
	Whether the refit is done on an expanding window including all the previous 
	data or a moving window where all previous data is used for the first estimation 
	and then moved by a length equal to refit.every (unless the window.size option 
	is used instead).}
  	\item{window.size}{
	If not NULL, determines the size of the moving window in the rolling estimation,
	which also determines the first point used.}
  	\item{solver}{ The solver to use.}
  	\item{fit.control}{ Control parameters parameters passed to the fitting 
  	function.}
  	\item{solver.control}{ Control parameters passed to the solver.}
  	\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
  	the parallel package. If it is not NULL, then this will be used for parallel 
  	estimation of the refits (remember to stop the cluster on completion).}
    \item{save.fit}{ Whether to save the fitted objects of class 
    \code{\linkS4class{DCCfit}} during the estimation of each (\dQuote{refit.every}). 
    If true, the directory to save must be provided. The function will not save 
    this by default for reasons of memory management, but can save it as an 
    \dQuote{.rda} file in the user's chosen directory for further analysis.}
  	\item{save.wdir}{ If \dQuote{save.fit} is true, the directory in which to 
  	save the \code{\linkS4class{DCCfit}} objects (1 for each \dQuote{refit.every}).}
  	\item{realizedVol}{ Required xts matrix for the realGARCH model.}
  	\item{\dots}{.}
}
\value{
A \code{\linkS4class{DCCroll}} object containing details of the DCC-GARCH 
rolling forecast.
}
\author{Alexios Ghalanos}
\keyword{methods}