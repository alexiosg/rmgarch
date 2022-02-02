\name{gogarchfit-methods}
\docType{methods}
\alias{gogarchfit}
\alias{gogarchfit,ANY-method}
\alias{gogarchfit,goGARCHspec-method}
\title{function: GO-GARCH Filter}
\description{
Method for filtering the GO-GARCH model.
}
\usage{
gogarchfit(spec, data, out.sample = 0, solver = "solnp", 
fit.control = list(stationarity = 1), solver.control = list(), cluster = NULL, 
VAR.fit = NULL, ARcoef = NULL, ...) 
}
\arguments{
\item{spec}{
A GO-GARCH spec object of class \code{\linkS4class{goGARCHspec}}.}
\item{data}{
A multivariate data object. Can be a matrix or data.frame or timeSeries.}
\item{out.sample}{
A positive integer indicating the number of periods before the last to keep for 
out of sample forecasting.}
\item{solver}{
One of either \dQuote{nlminb}, \dQuote{solnp} or \dQuote{gosolnp}. }
\item{solver.control}{
Control arguments list passed to optimizer.}
\item{fit.control}{
Control arguments passed to the fitting routine. Stationarity explicitly imposes 
the variance stationarity constraint during optimization.}
\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
the parallel package. If it is not NULL, then this will be used for parallel 
estimation (remember to stop the cluster on completion).}
\item{VAR.fit}{
(optional) A previously estimated VAR list returned from calling the 
\code{\link{varxfilter}} function.}
\item{ARcoef}{
An optional named matrix of the fitted AR parameters obtained from calling the 
\code{arfimafit} function on each series and then extracting the coefficients 
(the normal distribution should be used for the AR estimation). The number of  
columns should be equal to the number of series, and the rows should include the 
AR coefficients (common lag for all series), \sQuote{sigma}, and if included the 
mean (\sQuote{mu}). The option to pass the coefficients directly rather than 
letting the function estimate them may be useful for example when there are
convergence problems in the arfima routine and user control of each series 
estimation is desirable.
}
\item{...}{
Additional arguments passed to the ICA functions.}
}
\value{
A \code{\linkS4class{goGARCHfit}} object containing details of the GO-GARCH fit.
}
\note{
There is no check on the VAR.fit list passed to the method so particular care 
should be exercised so that the same data used in the fitting routine is also 
used in the VAR filter routine. The ability to pass this list of the 
pre-calculated VAR model is particularly useful when comparing different models 
(such as copula GARCH, DCC GARCH etc) using the same dataset and VAR method. 
Though the classical VAR estimation is very fast and may not require this extra 
step, the robust method is slow and therefore benefits from calculating this 
only once.
}
\author{Alexios Galanos}
\examples{
\dontrun{
data(dji30ret)
spec = gogarchspec(mean.model = list(demean = "constant"), 
variance.model = list(model = "sGARCH", garchOrder = c(1,1), submodel = NULL), 
distribution.model = list(distribution = "manig"),ica = "fastica")

fit = gogarchfit(spec = spec, data  = dji30ret[,1:4, drop = FALSE], 
out.sample = 50, gfun = "tanh")
fit
}}
\keyword{methods}
