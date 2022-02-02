\name{DCCtest}
\alias{DCCtest}
\title{
Engle and Sheppard Test of Dynamic Correlation
}
\description{
A test of non-constant correlation based on Engle and Sheppard (2001).
}
\usage{
DCCtest(Data, garchOrder = c(1,1), n.lags = 1, solver = "solnp", 
solver.control = list(), cluster = NULL, Z = NULL)
}
\arguments{
\item{Data}{
A multivariate data matrix.
}
\item{garchOrder}{
The first stage common GARCH order.
}
\item{n.lags}{
The number of lags to test for the presence of non-constant correlation.
}
\item{solver}{ Either \dQuote{solnp} or \dQuote{nlminb} .}
\item{solver.control}{Control arguments list passed to optimizer.}
\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
the parallel package. If it is not NULL, then this will be used for parallel 
estimation (remember to stop the cluster on completion).}
\item{Z}{
(Optional) The standardized residuals from a constant correlation model. If 
supplied the model is not estimated since this is the only input the test
requires.
}
}
\details{
The test effectively equates to estimating a multivariate dataset using the
Constant Conditional Correlation (CCC) model of Bollerslev (1990) and after 
which the standardized residuals (standardized by the symmetric square root 
decomposition of the estimated constant correlation matrix) should be i.i.d. 
with covariance the identity matrix. Testing for this can be done using a 
series of artificial regressions on the outer and lagged product of these 
residuals and a constant. In the rmgarch package, the CCC model is calculated
using a static GARCH copula (Normal) model.
}
\value{
A list with the proposed Null hypothesis (H0), the test statistic and its 
p-value.
}
\author{Alexios Galanos}
\references{
Bollerslev, T. 1990, Modelling the coherence in short-run nominal exchange 
rates: a multivariate generalized ARCH model, \emph{The Review of Economics and 
Statistics}, \bold{72(3)}, 498--505.\cr
Engle, R.F. and Sheppard, K. 2001, Theoretical and empirical properties of 
dynamic conditional correlation multivariate GARCH, \emph{NBER Working Paper}.\cr
}
