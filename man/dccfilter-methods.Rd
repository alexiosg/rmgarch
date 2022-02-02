\name{dccfilter-methods}
\docType{methods}
\alias{dccfilter}
\alias{dccfilter,ANY-method}
\alias{dccfilter,DCCspec-method}
\title{function: DCC-GARCH Filter}
\description{
Method for creating a DCC-GARCH filter object.
}
\usage{
dccfilter(spec, data, out.sample = 0, filter.control = list(n.old = NULL), 
cluster = NULL, varcoef = NULL, realizedVol = NULL, ...) 
}
\arguments{
  	\item{spec}{ A \code{\linkS4class{DCCspec}} object created by calling 
  	\code{\link{dccspec}} with fixed parameters for the coefficients.}
  	\item{data}{ A multivariate data object of class xts, or one which can be
  	coerced to such.}
  	\item{out.sample}{ A positive integer indicating the number of periods 
  	before the last to keep for out of sample forecasting.}
	\item{filter.control}{Control arguments passed to the filtering routine (see
	note.}
	\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
  	the parallel package. If it is not NULL, then this will be used for parallel 
  	estimation (remember to stop the cluster on completion).}
	\item{varcoef}{ If a VAR model was chosen, then this is the VAR coefficient 
	matrix which must be supplied. No checks are done on its dimension or 
	correctness so it is up to the user to perform the appropriate checks.}
	\item{realizedVol}{ Required xts matrix for the realGARCH model.}
	\item{...}{ . }
}
\value{
A \code{\linkS4class{DCCfilter}} object containing details of the DCC-GARCH 
filter.
}
\note{
The \sQuote{n.old} option in the \code{filter.control} argument is key in 
replicating conditions of the original fit. That is, if you want to filter a 
dataset consisting of an expanded dataset (versus the original used in fitting), 
but want to use the same assumptions as the original dataset then the \sQuote{n.old} 
argument denoting the original number of data points passed to the 
\code{\link{dccfit}} function must be provided. This is then used to ensure 
that some calculations which make use of the full dataset (unconditional 
starting values for the garch filtering and the dcc model) only use the first 
\sQuote{n.old} points thus replicating the original conditions making filtering 
appropriate for rolling 1-ahead forecasting.\cr
For extensive examples look in the \sQuote{rmgarch.tests} folder.
}
\author{Alexios Galanos}
\keyword{methods}
