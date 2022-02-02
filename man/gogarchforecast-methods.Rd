\name{gogarchforecast-methods}
\docType{methods}
\alias{gogarchforecast}
\alias{gogarchforecast-methods}
\alias{gogarchforecast,ANY-method}
\alias{gogarchforecast,goGARCHfit-method}
\title{function: GO-GARCH Forecast}
\description{
Method for forecasting from the GO-GARCH model.
}
\usage{
gogarchforecast(fit, n.ahead = 10, n.roll = 0, 
external.forecasts = list(mregfor = NULL), cluster = NULL, ...) 
}
\arguments{
\item{fit}{
A GO-GARCH fit object of class \code{\linkS4class{goGARCHfit}}.}
\item{n.ahead}{
The forecast horizon.}
\item{n.roll}{
The no. of rolling forecasts to create beyond the first one.}
\item{external.forecasts}{
A list with a matrix object of the external lagged forecasts (if used). These 
must contain (n.roll+1) x n.ahead forecasts. }
\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
the parallel package. If it is not NULL, then this will be used for parallel 
estimation (remember to stop the cluster on completion).}
\item{...}{.}
}
\value{
A \code{\linkS4class{goGARCHforecast}} object containing details of the GO-GARCH 
forecast.
}
\author{Alexios Galanos}
\examples{
\dontrun{
data(dji30ret)
spec = gogarchspec()
fit = gogarchfit(spec = spec, data = dji30ret[,1:4], out.sample = 10, 
gfun = "tanh")
forecast = gogarchforecast(fit, n.ahead = 1, n.roll = 9)
}}
\keyword{methods}
