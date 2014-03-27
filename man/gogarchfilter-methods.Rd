\name{gogarchfilter-methods}
\docType{methods}
\alias{gogarchfilter}
\alias{gogarchfilter,ANY-method}
\alias{gogarchfilter,goGARCHfit-method}
\title{function: GO-GARCH Filter}
\description{
Method for filtering the GO-GARCH model.
}
\usage{
gogarchfilter(fit, data, out.sample = 0, n.old = NULL, cluster = NULL, ...) 
}
\arguments{
\item{fit}{
A GO-GARCH fit object of class \code{\linkS4class{goGARCHfit}}.}
\item{data}{
A multivariate data object. Can be a matrix or data.frame or timeSeries.}
\item{out.sample}{
A positive integer indicating the number of periods before the last to keep for 
out of sample forecasting.}
\item{n.old}{
For comparison with goGARCHfit models using the out.sample argument, this is the 
length of the original dataset.}
\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
the parallel package. If it is not NULL, then this will be used for parallel 
estimation (remember to stop the cluster on completion).}
\item{...}{.}
}
\value{
A \code{\linkS4class{goGARCHfilter}} object containing details of the GO-GARCH 
filter.
}
\author{Alexios Ghalanos}
\examples{
\dontrun{
data(dji30ret)
spec = gogarchspec()
fit = gogarchfit(spec = spec, data = dji30ret[,1:4], gfun = "tanh")
filter = gogarchfilter(fit, data = dji30ret[,1:4])
}}
\keyword{methods}