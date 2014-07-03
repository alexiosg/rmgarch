\name{cgarchfilter-methods}
\docType{methods}
\alias{cgarchfilter}
\alias{cgarchfilter,ANY-method}
\alias{cgarchfilter,cGARCHspec-method}
\alias{cgarchfilter,cGARCHfit-method}
\alias{cgarchfilter,cGARCHfilter-method}

\title{function: Copula-GARCH Filter}
\description{
Method for creating a Copula-GARCH filter object.
}
\usage{
cgarchfilter(mod, data, newdata = NULL, out.sample = 0, 
filter.control = list(n.old = NULL), spd.control = list(lower = 0.1, 
upper = 0.9, type = "pwm", kernel = "epanech"), 
cluster = NULL, varcoef = NULL, realizedVol = NULL, ...)  
}
\arguments{
  	\item{mod}{ Either a \code{\linkS4class{cGARCHspec}} object created by 
  	calling \code{\link{cgarchspec}} with fixed parameters for the coeffficients,
  	a \code{\linkS4class{cGARCHfit}} or \code{\linkS4class{cGARCHfilter}} object.}
  	\item{data}{ A multivariate xts data object or one which can be coerced to
  	such.}
  	\item{newdata}{ In the case of 	a \code{\linkS4class{cGARCHfit}} or 
  	\code{\linkS4class{cGARCHfilter}} object, instead of filtering based on the
  	entire dataset (data), the user can pass only the data which is new which is
  	then filtered based on the last values of the existing dataset in the object.
  	The returned object will then append this new data to the old data so that
  	the returned object can then be passed again to this method with updated data.
  	This provides the possibility of filtering in almost real-time.}
  	\item{out.sample}{ A positive integer indicating the number of periods 
  	before the last to keep for out of sample forecasting.}
	\item{filter.control}{Control arguments passed to the filtering routine (see
	note below).}
  	\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
  	the parallel package. If it is not NULL, then this will be used for parallel 
  	estimation (remember to stop the cluster on completion).}
    \item{spd.control}{ If the spd transformation was chosen in the 
	specification, the spd.control passes its arguments to the 
	\code{spdfit} routine of the \code{spd} package.}
	\item{varcoef}{ If a VAR model was chosen, then this is the VAR coefficient 
	matrix which must be supplied. No checks are done on its dimension or 
	correctness so it is up to the user to perform the appropriate checks.}
	\item{realizedVol}{ Required xts matrix for the realGARCH model.}
	\item{...}{ . }
}
\value{
A \code{\linkS4class{cGARCHfilter}} object containing details of the 
Copula-GARCH filter and sharing most of the methods of the 
\code{\linkS4class{cGARCHfit}} class.
}
\note{
The \sQuote{n.old} option in the \code{filter.control} argument is key in 
replicating conditions of the original fit. That is, if you want to filter a 
dataset consisting of an expanded dataset (versus the original used in fitting), 
but want to use the same assumptions as the original dataset then the \sQuote{n.old} 
argument denoting the original number of data points passed to the 
\code{\link{cgarchfit}} function must be provided. This is then used to ensure 
that some calculations which make use of the full dataset (unconditional 
starting values for the garch filtering, the dcc model and the copula 
transformation methods) only use the first \sQuote{n.old} points thus 
replicating the original conditions making filtering appropriate for rolling 
1-ahead forecasting.\cr
For extensive examples look in the \sQuote{rmgarch.tests} folder.
}
\author{Alexios Ghalanos}
\keyword{methods}