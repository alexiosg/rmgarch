\name{fmoments-methods}
\docType{methods}
\alias{fmoments}
\alias{fmoments-methods}
\alias{fmoments,ANY-method}
\title{Moment Based Forecast Generation}
\description{
Generates n-ahead forecast moment matrices given a choice of data generating 
processes.
}
\usage{
fmoments(spec, Data, n.ahead = 1, roll  = 0, solver = "solnp", 
solver.control = list(), fit.control = list(eval.se = FALSE), 
cluster = NULL, save.output = FALSE, save.dir = getwd(), 
save.name = paste("M", sample(1:1000, 1), sep = ""), ...)
}
\arguments{
\item{Data}{ An n-by-m data matrix or data.frame.}
\item{spec}{ Either a DCCspec or GOGARCHspec.}
\item{n.ahead}{ The n.ahead forecasts (n.ahead>1 is unconditional).}
\item{roll}{ Whether to fit the data using (n - roll) periods and then
return a (roll+1) n-ahead rolling forecast moments.}

\item{solver}{ The choice of solver to use for all models but \dQuote{var},
and includes \sQuote{solnp}, \sQuote{nlminb} and \sQuote{nloptr}.}
\item{solver.control}{ Optional control options passed to the appropriate solver
chosen.}
\item{fit.control}{Control arguments passed to the fitting routine.}
\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
the parallel package. If it is not NULL, then this will be used for parallel 
estimation of the refits (remember to stop the cluster on completion).}
\item{save.output}{ Whether output should be saved to file instead of being 
returned to the workspace.}
\item{save.dir}{ The directory to save output if save.output is TRUE.}
\item{save.name}{ The name of the file to save the output list.}
\item{...}{ Additional parameters passed to the model fitting routines. In 
particular, for the \sQuote{gogarch} model additional parameters are passed to
the ICA routines, whereas for the \sQuote{dcc} and \sQuote{cgarch} models 
this would include the \sQuote{realizedVol} xts matrix for the realGARCH model.}
}
\value{
A \code{\linkS4class{fMoments}} object containing the forecast moments 
(list of length roll+1) and the model details (list).
}
\details{
The function allows to generate forecast covariance matrices for use in the
QP based EV model, and also for the \dQuote{gogarch} model higher co-moment 
matrices for use in the Utility maximization model implemented separately.
}
\author{Alexios Galanos}
\keyword{methods}
