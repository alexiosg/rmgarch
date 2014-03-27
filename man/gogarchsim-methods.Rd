\name{gogarchsim-methods}
\docType{methods}
\alias{gogarchsim}
\alias{gogarchsim,ANY-method}
\alias{gogarchsim,goGARCHfit-method}
\title{function: GO-GARCH Simulation}
\description{
Method for simulation from a fitted GO-GARCH model.
}
\usage{
gogarchsim(fit, n.sim = 1, n.start = 0, m.sim = 1, 
startMethod = c("unconditional", "sample"), prereturns = NA, preresiduals = NA,
presigma = NA, mexsimdata = NULL, rseed = NULL, cluster = NULL, ...)
}
\arguments{
\item{fit}{
A GO-GARCH fit object of class \code{\linkS4class{goGARCHfit}}.}
\item{n.sim}{
The simulation horizon.}
\item{n.start}{
The burn-in sample.}
\item{m.sim}{
The number of simulations.}
\item{startMethod}{
Starting values for the simulation. Valid methods are \dQuote{unconditional} 
for the  expected values given the density, and \dQuote{sample} for the 
ending values of the actual data from the fit object.}
\item{prereturns}{
Allows the starting return data to be provided by the user.}
\item{preresiduals}{
Allows the starting factor residuals to be provided by the user.}
\item{presigma}{
Allows the starting conditional factor sigma to be provided by the user.}
\item{mexsimdata}{
A list of matrices with the simulated lagged external variables (if any). The 
list should be of size m.sim and the matrices each have n.sim + n.start rows.}
\item{rseed}{
Optional seeding value(s) for the random number generator.}
\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
the parallel package. If it is not NULL, then this will be used for parallel 
estimation (remember to stop the cluster on completion).}
\item{\dots}{.}
}
\value{
A \code{\linkS4class{goGARCHsim}} object containing details of the GO-GARCH 
simulation.
}
\author{Alexios Ghalanos}
\keyword{methods}