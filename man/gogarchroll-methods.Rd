\name{gogarchroll-methods}
\alias{gogarchroll}
\alias{gogarchroll,ANY-method}
\alias{gogarchroll,goGARCHspec-method}
\title{function: GO-GARCH Rolling Estimation}
\description{
Method for performing rolling estimation of the GO-GARCH model.
}
\usage{
gogarchroll(spec, data, n.ahead = 1, forecast.length = 50, n.start = NULL, 
refit.every = 25, refit.window = c("recursive", "moving"), window.size = NULL, 
solver = "solnp", solver.control = list(), fit.control = list(), rseed = NULL,  
cluster = NULL, save.fit = FALSE, save.wdir = NULL, ...)
}
\arguments{
\item{spec}{
A GO-GARCH spec object of class \code{\linkS4class{goGARCHspec}}.}
\item{data}{
A multivariate data object. Can be a matrix or data.frame or timeSeries. }
\item{n.ahead}{
The forecast horizon (only 1-ahead supported for rolling forecasts). }
\item{forecast.length}{
The length of the total forecast for which out of sample data from the dataset 
will be excluded for testing. }
\item{n.start}{ Instead of forecast.length, this determines the starting 
point in the dataset from which to initialize the rolling forecast.}
\item{refit.every}{
Determines every how many periods the model is re-estimated. }
\item{refit.window}{
Whether the refit is done on an expanding window including all the previous data 
or a moving window where all previous data is used for the first estimation 
and then moved by a length equal to refit.every (unless the window.size option 
is used instead).
}
\item{window.size}{
If not NULL, determines the size of the moving window in the rolling estimation,
which also determines the first point used.
}
\item{solver}{
The solver to use. }
\item{fit.control}{
Control parameters parameters passed to the fitting function. }
\item{solver.control}{
Control parameters passed to the solver. }
\item{rseed}{
Initialization seed for first ICA fit. The rest of the ICA fits are initialized 
with the previous mixing matrix (using A.init). }
\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
the parallel package. If it is not NULL, then this will be used for parallel 
estimation (remember to stop the cluster on completion).}
\item{save.fit}{
Whether to save the fitted objects of class \code{\linkS4class{goGARCHfit}} 
during the estimation of each (\dQuote{refit.every}). If true, the directory to 
save must be provided (see below). The function will not save this by default 
for reasons of memory management, but can save it as an \dQuote{.rda} file in 
the user's chosen directory for further analysis.}
\item{save.wdir}{
If \dQuote{save.fit} is true, the directory in which to save the 
\code{\linkS4class{goGARCHfit}} objects (1 for each \dQuote{refit.every}).}
\item{\dots}{.}
}
\value{
An object of class \code{\linkS4class{goGARCHroll}}.
}
\author{Alexios Ghalanos}
\keyword{methods}