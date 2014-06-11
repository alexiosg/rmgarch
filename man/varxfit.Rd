\name{varxfit}
\alias{varxfilter}
\alias{varxfit}
\alias{varxforecast}
\alias{varxsim}
\title{
VARX Fit/Filter/Forecast/Simulation Functions
}
\description{
Vector Autoregressive (VAR) with Constant and Optional Exogenous Regressors (X)
Fit, Filter, Forecast and Simulation functions for use with multivariate GARCH 
models.
}
\usage{
varxfit(X, p, constant = TRUE, exogen = NULL, robust = FALSE, gamma = 0.25, delta = 0.01, nc = 10, 
ns = 500, postpad = c("none", "constant", "zero", "NA"), cluster = NULL)
varxfilter(X, p, Bcoef, exogen = NULL, postpad = c("none", "constant", "zero", "NA"))
varxforecast(X, Bcoef, p, out.sample, n.ahead, n.roll, mregfor)
varxsim (X, Bcoef, p, n.sim, n.start, prereturns, resids, mexsimdata)
}
\arguments{
  \item{X}{
  A multivariate data matrix.
}
  \item{p}{
  The number of autoregressive lags.
}
  \item{constant}{
	Whether to include a constant.
}
  \item{exogen}{
  An optional matrix of exogenous regressors with as many rows as X, and 
  appropriately lagged.
}
  \item{Bcoef}{
  A matrix of coefficients for the varxfilter function.
}
  \item{robust}{
  Whether to use the robust version of VAR based on the multivariate Least 
  Trimmed Squares Estimator described in Croux and Joossens (2008).
}
 \item{gamma}{
 Proportion to trim in the robust method.
 }
  \item{delta}{
 The critical value for Reweighted estimator for the robust method.
 }
  \item{ns}{
 The number of subsets to use for the robust method.
 }
  \item{nc}{
 The number of C-steps to use for the robust method.
 }
 \item{postpad}{
 (defaults to \sQuote{none}) Whether to postpad the fitted/filtered values 
 (and hence calculation of residuals) with the estimated constant, zeros or NA's, 
 thus returning matrices of the same size as the input data (rather than input 
 data size less the number of lags).
 }
 \item{cluster}{ A cluster object created by calling \code{makeCluster} from 
 the parallel package. If it is not NULL, then this will be used for parallel 
 estimation in the case of the robust VAR version (remember to stop the cluster 
 on completion).}
 \item{out.sample}{
 The number of points kept for out of sample rolling forecast.
 }
 \item{n.ahead}{
 The forecast horizon.
 }
 \item{n.roll}{
 In combination with \code{out.sample}, determines the number of times to roll
 forward the n.ahead forecast using data left out of sample.
 }
 \item{mregfor}{
 Matrix of external regressor forecasts (with appropriate lag structure).
 }
 \item{n.sim}{
 Simulation horizon.
 }
 \item{n.start}{
 Simulation burn-in sample.
 }
 \item{prereturns}{
 Optionally supplied pre-return matrix with \dQuote{p} lags to initialize 
 simulation.
 }
 \item{resids}{
 Matrix of randomly generated residuals of size n.sim+n.start.
 }
 \item{mexsimdata}{
 Matrix of external regressor pre-generated random values to use in the
 simulation (if NULL then assumed zero).
 }
}
\details{
This are convenience functions to be optionally used when using the multivariate 
GARCH methods.
}
\value{
A list with the following items:\cr
\item{Bcoef }{[varxfit, varxfilter]\cr
The coefficient matrix with rows equal to number of assets, and columns equal to
number of assets x number of lags plus 1 (constant) plus number of exogenous 
regressors.}
\item{xfitted }{[varxfit, varxfilter]\cr
The fitted/filtered series (conditional mean series).}
\item{xresiduals }{[varxfit, varxfilter]\cr
The residuals.}
\item{Bcov }{[varxfit]\cr
The covariance matrix of the coefficients.}
\item{se }{[varxfit]\cr
The standard error of the coefficients.} 
\item{tstat }{[varxfit]\cr
The t-stat of the s.e.}
\item{pstat }{[varxfit]\cr
The p-values of the s.e.}
\item{lag }{[varxfit, varxfilter]\cr
The number of autoregressive lags.}
\item{mxn}{[varxfit]\cr
The number of exogenous regressors .}
\item{meansim}{[varxsim]\cr
The simulated conditional mean.}\cr
The varxsim returns an n.sim x n.asset matrix of the simulated conditional
means, while the varxforecast returns an n.ahead x n.assets x (n.roll+1) array 
of the forecast conditional means.
}
\references{
Lutkepohl, H. 2005, New introduction to multiple time series analysis, 
\emph{Springer}.\cr
Croux, C. and Joossens, K. 2008, Robust estimation of the vector autoregressive 
model by a least trimmed squares procedure, \emph{COMPSTAT}, 489--501.\cr
}
\author{Alexios Ghalanos}
\note{
Part of the varxfit functionality and structure is inspired from the 'vars' 
package, but the estimation method is implemented in a very quick way without 
calling 'lm'. The robust method is based on the matlab program of Christophe 
Croux available from his website and the option of using parallel computation is
implemented for this particular choice.The \code{postpad} option is used when the
returned data needs to be of the same size as the inputed data for easier
manipulation/comparison (since padding is done post-estimation, there is no bias 
introduced during estimation).
}