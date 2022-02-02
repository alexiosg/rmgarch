\name{cgarchspec-methods}
\docType{methods}
\alias{cgarchspec}
\alias{cgarchspec,ANY-method}
\alias{cgarchspec,uGARCHmultispec-method}
\title{function: Copula-GARCH Specification}
\description{
Method for creating a Copula-GARCH specification object prior to fitting.
}
\usage{
cgarchspec(uspec, VAR = FALSE, robust = FALSE, lag = 1, lag.max = NULL, 
lag.criterion = c("AIC", "HQ", "SC", "FPE"), external.regressors = NULL, 
robust.control = list(gamma = 0.25, delta = 0.01, nc = 10, ns = 500), 
dccOrder = c(1, 1), asymmetric = FALSE, 
distribution.model = list(copula = c("mvnorm", "mvt"), 
method = c("Kendall", "ML"), time.varying = FALSE, 
transformation = c("parametric", "empirical", "spd")), 
start.pars = list(), fixed.pars = list()) 
}
\arguments{
  \item{uspec}{
  A \code{\linkS4class{uGARCHmultispec}} object created by calling 
  \code{\link{multispec}} on a list of univariate GARCH specifications.}
  \item{VAR}{ Whether to fit a VAR model for the conditional mean.}
  \item{robust}{ Whether to use the robust version of VAR.}
  \item{lag}{ The VAR lag.}
  \item{lag.max}{ The maximum VAR lag to search for best fit.}
  \item{lag.criterion}{ The criterion to use for choosing the best lag when
  lag.max is not NULL.}
  \item{external.regressors}{ Allows for a matrix of common pre-lagged external 
  regressors for the VAR option.
  }
  \item{robust.control}{ The tuning parameters to the robust regression 
  including the proportion to trim (\dQuote{gamma}), the critical value for 
  reweighted estimator (\dQuote{delta}), the number of subsets (\dQuote{ns}) and 
  the number of C-steps (\dQuote{nc}.}
  \item{dccOrder}{ The DCC autoregressive order.}
  \item{asymmetric}{ Whether to include an asymmetry term to the DCC model (thus
  estimating the aDCC).}
  \item{distribution.model}{ The Copula distribution model. Currently the 
  multivariate Normal and Student Copula are supported. }
  \item{time.varying}{ Whether to fit a dynamic DCC Copula. }
  \item{transformation}{ The type of transformation to apply to the marginal 
  innovations of the GARCH fitted models. Supported methods are parametric 
  (Inference Function of Margins), empirical (Pseudo ML), and Semi-Parametric 
  using a kernel interior and GPD tails (via the \sQuote{spd} package).}
  \item{start.pars}{ (optional) Starting values for the DCC parameters 
  (starting values for the univariate garch specification should be passed 
  directly via the \sQuote{uspec} object).}
  \item{fixed.pars}{ (optional) Fixed DCC parameters.}
}
\value{
A \code{\linkS4class{cGARCHspec}} object containing details of the Copula-GARCH 
specification.
}
\details{
The transformation method allows for parametric (Inference-Functions for Margins), 
empirical (Pseudo-Likelihood) and semi-parametric (via the spd package).\cr
When the Student Copula is jointly estimated with student margins having so that 
a common shape parameter is obtained, this results in the multivariate Student 
distribution. When estimating the Student Copula with disparate margins, a 
meta-student distribution is obtained. Additionally, the correlation parameter 
in the static Student Copula may be estimated either by Kendall's tau 
transformation or Maximum Likelihood.\cr
The \code{robust} option allows for a robust version of VAR based on the 
multivariate Least Trimmed Squares Estimator described in Croux and Joossens 
(2008).
}
\author{Alexios Galanos}
\keyword{methods}
