\name{dccspec-methods}
\docType{methods}
\alias{dccspec}
\alias{dccspec,ANY-method}
\alias{dccspec,uGARCHmultispec-method}
\title{function: DCC-GARCH Specification}
\description{
Method for creating a DCC-GARCH specification object prior to fitting.
}
\usage{
dccspec(uspec, VAR = FALSE, robust = FALSE, lag = 1, lag.max = NULL,
lag.criterion = c("AIC", "HQ", "SC", "FPE"), external.regressors = NULL,
robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500),
dccOrder = c(1,1), model = c("DCC", "aDCC", "FDCC"), groups = rep(1, length(uspec@spec)),
distribution = c("mvnorm", "mvt", "mvlaplace"), start.pars = list(), fixed.pars = list())
}
\arguments{
  \item{uspec}{
  A \code{uGARCHmultispec} object from the \pkg{rugarch} package, created by calling
  \code{\link[rugarch]{multispec}} on a list of univariate GARCH specifications.}
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
  re-weighted estimator (\dQuote{delta}), the number of subsets (\dQuote{ns})
  and the number of C-steps (\dQuote{nc}.}
  \item{dccOrder}{ The DCC autoregressive order.}
  \item{model}{ The DCC model to use, with a choice of the symmetric DCC,
  asymmetric (aDCC) and the Flexible DCC (FDCC). See notes for more details.}
  \item{groups}{ The groups corresponding to each asset in the FDCC model, where
  these are assumed and checked to be contiguous and increasing (unless only 1 group).}
  \item{distribution}{ The multivariate distribution. Currently the multivariate
  Normal, Student and Laplace are implemented, and only the Normal for the FDCC model.}
  \item{start.pars}{ (optional) Starting values for the DCC parameters (starting
  values for the univariate garch specification should be passed directly
  via the \sQuote{uspec} object).}
  \item{fixed.pars}{ (optional) Fixed DCC parameters. This is required in the
  \code{\link{dccfilter}}, \code{\link{dccforecast}}, \code{\link{dccsim}} with
  spec, and  \code{\link{dccroll}}
  methods.}
}
\value{
A \code{\linkS4class{DCCspec}} object containing details of the DCC-GARCH
specification.
}
\details{
The \code{robust} option allows for a robust version of VAR based on the
multivariate Least Trimmed Squares Estimator described in Croux and Joossens
(2008).
}
\note{
The FDCC model of Billio, Caporin and Gobbo (2006) allows different DCC
parameters to govern the dynamics of the correlation of distinct groups. The
drawback is a somewhat larger parameter set, and no correlation targeting.
Still, it remains a feasible model for not too large a number of groups, and
avoids the unrealistic assumption, particularly for large datasets, of one
parameter governing all the dynamics, as in the DCC model. Note that the group
indices must be increasing (unless all 1), which means that you should arrange
your dataset so that the assets are ordered by their groups.
}
\references{
Billio, M., Caporin, M., & Gobbo, M. 2006, Flexible dynamic conditional
correlation multivariate GARCH models for asset allocation, \emph{Applied
Financial Economics Letters}, \bold{2(02)}, 123--130.\cr
Croux, C. and Joossens, K. 2008, Robust estimation of the vector autoregressive
model by a least trimmed squares procedure, \emph{COMPSTAT}, 489--501.\cr
Cappiello, L., Engle, R.F. and Sheppard, K. 2006, Asymmetric dynamics in the
correlations of global equity and bond returns, \emph{Journal of Financial
Econometrics} \bold{4}, 537--572.\cr
Engle, R.F. and Sheppard, K. 2001, Theoretical and empirical properties of
dynamic conditional correlation multivariate GARCH, \emph{NBER Working Paper}.\cr
}
\author{Alexios Galanos}
\keyword{methods}
