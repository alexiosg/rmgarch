\name{goGARCHfit-class}
\docType{class}
\alias{goGARCHfit-class}
\alias{as.matrix,goGARCHfit-method}
\alias{coef,goGARCHfit-method}
\alias{convolution,goGARCHfit-method}
\alias{convolution}
\alias{fitted,goGARCHfit-method}
\alias{residuals,goGARCHfit-method}
\alias{likelihood,goGARCHfit-method}
\alias{nisurface,goGARCHfit-method}
\alias{nisurface}
\alias{gportmoments,goGARCHfit-method}
\alias{gportmoments}
\alias{rcor,goGARCHfit-method}
\alias{rcor}
\alias{rcoskew,goGARCHfit-method}
\alias{rcoskew}
\alias{rcokurt,goGARCHfit-method}
\alias{rcokurt}
\alias{rcov,goGARCHfit-method}
\alias{rcov}
\alias{betacovar}
\alias{betacovar,goGARCHfit-method}
\alias{betacoskew}
\alias{betacoskew,goGARCHfit-method}
\alias{betacokurt}
\alias{betacokurt,goGARCHfit-method}
\alias{show,goGARCHfit-method}
\title{class: GO-GARCH Fit Class}
\description{Class for the GO-GARCH fitted object.}
\section{Objects from the Class}{
The class is returned by calling the function \code{\link{gogarchfit}}.
}
\section{Slots}{
\describe{
\item{\code{mfit}:}{Multivariate fit object.}
\item{\code{model}:}{Object of class \code{"vector"} containing details of the
GO-GARCH model specification.}
}
}
\section{Extends}{
Class \code{mGARCHfit}, directly.\\
Class \code{GARCHfit} object from the \pkg{rugarch} package, by class \code{mGARCHfit}, distance 2.\\
Class \code{rGARCH} object from the \pkg{rugarch} package, by class \code{mGARCHfit}, distance 3.
}
\section{Methods}{
\describe{
\item{as.matrix}{\code{signature(x = "goGARCHfit")}:
function:\cr
\bold{as.matrix(x, which = "A")}\cr
This returns four types of matrices relating to the estimation of the
independent components in the GO-GARCH model. Valid choices are \dQuote{A}
for the mixing matrix, \dQuote{W} for the unmixing matrix, \dQuote{U} for the
rotational matrix and \dQuote{K} for the whitening matrix, \dQuote{Kinv} for
the de-whitening matrix.}
\item{coef}{\code{signature(object = "goGARCHfit")}:  extraction of independent
factor GARCH model coefficients.}
\item{convolution}{\code{signature(object = "goGARCHfit")}:\cr
function:\cr
\bold{convolution(object, weights, fft.step = 0.001, fft.by = 0.0001,
fft.support = c(-1, 1), support.method = c("user", "adaptive"), use.ff = TRUE,
cluster = NULL, trace = 0,...)}\cr
The convolution method takes a goGARCHfit object and a weights vector or
matrix and calculates the weighted density. If a vector is given, it must be
the same length as the number of assets, otherwise a matrix with row
dimension equal to the row dimension of the filtered dataset (i.e. less any
lags). In the case of the multivariate normal distribution, this simply
returns the linear and quadratic transformation of the mean and covariance
matrix, while in the multivariate affine NIG distribution this is based on
the numerical inversion by FFT of the characteristic function. In that case,
the \dQuote{fft.step} option determines the stepsize for tuning the
characteristic function inversion, \dQuote{fft.by} determines the resolution
for the equally spaced support given by \dQuote{fft.support}, while the use
of the \dQuote{ff} package is recommended to avoid memory problems on some
systems and is turned on via the \dQuote{use.ff} option.
The \dQuote{support.method} option allows either a fixed support range to be
given (option \sQuote{user}), else an adaptive method is used based on the
min and max of the assets at each point in time at the 0.00001 and 1-0.00001
quantiles. The range is equally spaced subject to the \dQuote{fft.by} value
but the returned object no longer makes of the \dQuote{ff} package returning
instead a list. Finally, the option for parallel computation is available via
the use of a cluster object as elsewhere in this package.}
\item{fitted}{\code{signature(object = "goGARCHfit")}:
Extracts the conditional mean equation fitted values.}
\item{residuals}{\code{signature(object = "goGARCHfit")}:
Extracts the conditional mean equation residuals.}
\item{likelihood}{\code{signature(object = "goGARCHfit")}:
The quasi log-likelihood of the model, which being an independent factor model
is the sum of the univariate GARCH log-likelihoods plus a term for the mixing
matrix. For a dimensionality reduced system, this is NA.}
\item{nisurface}{\code{signature(object = "goGARCHfit")}:\cr
function:\cr
\bold{nisurface(object, type = "cov", pair = c(1, 2), factor = c(1,2),
plot = TRUE)}\cr
Creates the covariance or correlation (determined by \dQuote{type} being
either \dQuote{cov} or \dQuote{cor}) news impact surface for a pair of
assets and factors. Since the shocks impact the factors independently,
the news impact surface is a combination of the independent news impact
curves of the factors which when combined via the mixing matrix A, create
the dynamics for the underlying asset-factor surface function.}
\item{gportmoments}{\code{signature(object = "goGARCHfit")}:\cr
function:\cr
\bold{gportmoments(object, weights)}\cr
Calculates the first 4 portfolio moments using the geometric properties of
the model, given a vector or matrix of asset weights. If a matrix is given
it must have row dimension equal to the row dimension of the filtered
dataset (i.e. less any lags), else if a vector is given it will be replicated
for all time points.}
\item{rcoskew}{\code{signature(object = "goGARCHfit")}
function:\cr
\bold{rcoskew(object, standardize = TRUE, from = 1, to = 1)}\cr
Returns the 'time-varying'  NxN^2 coskewness tensor in array format. The
\dQuote{from} and \dQuote{to} options indicate the time indices for which
to return the arrays. Because of memory issues, this is limited to 100
indices per call.}
\item{rcokurt}{\code{signature(object = "goGARCHfit")}
function:\cr
\bold{rcokurt(object, standardize = TRUE, from = 1, to = 1)}\cr
Returns the 'time-varying'  NxN^3 cokurtosis tensor in array format. The
\dQuote{from} and \dQuote{to} options indicate the time indices for which
to return the arrays. Because of memory issues, this is limited to models
with less than 100 assets.}
\item{rcov}{\code{signature(object = "goGARCHfit")}:
Returns the time-varying NxN covariance matrix in array format unless
\sQuote{output} is set to \dQuote{matrix} in which case the array
is flattened and the lower and main diagonal time varying values are
returned (and if a date exists, then the returned object is of class xts).
}
\item{rcor}{\code{signature(object = "goGARCHfit")}:
Returns the time-varying NxN correlation matrix in array format unless
\sQuote{output} is set to \dQuote{matrix} in which case the array
is flattened and the lower and main diagonal time varying values are
returned (and if a date exists, then the returned object is of class xts).
}
\item{betacovar}{\code{signature(object = "goGARCHfit")}:
function:\cr
\bold{betacovar(object, weights, asset = 1, cluster = NULL)}\cr
Returns the covariance beta given a matrix (of length equal to the number of
rows of the original data, or vector which is then recycled to the number
of rows of the original data) of benchmark weights and the asset number.}
\item{betacoskew}{\code{signature(object = "goGARCHfit")}:
function:\cr
\bold{betacoskew(object, weights, asset = 1, cluster = NULL)}\cr
Returns the coskewness beta given a matrix (of length equal to the number of
rows of the original data, or vector which is then recycled to the number
of rows of the original data) of benchmark weights and the asset number.}
\item{betacokurt}{\code{signature(object = "goGARCHfit")}:
function:\cr
\bold{betacokurt(object, weights, asset = 1, cluster = NULL)}\cr
Returns the cokurtosis beta given a matrix (of length equal to the number of
rows of the original data, or vector which is then recycled to the number
of rows of the original data) of benchmark weights and the asset number.}
\item{show}{\code{signature(object = "goGARCHfit")}: Summary method.}
}
}
\references{
de Athayde, G.M. and Flores Jr, R.G. 2002, On Certain Geometric Aspects of
Portfolio Optimisation with Higher Moments, \emph{mimeo}.\cr
Broda, S.A. and Paolella, M.S. 2009, CHICAGO: A Fast and Accurate Method for
Portfolio Risk Calculation, \emph{Journal of Financial Econometrics} \bold{7(4)},
412--436.\cr
Paolella, M.S. 2007, Intermediate Probability - A Computational Approach,
\emph{Wiley-Interscience}.\cr
Schmidt, R., Hrycej, T., and Stutzle (2006), Multivariate distribution models with
generalized hyperbolic margins, \emph{Computational Statistics & Data Analysis}
\bold{50(8)}, 2065---2096.\cr
}
\author{Alexios Galanos}
\note{
The reference by Paolella (2007) contains more details on the algorithm for the
characteristic function inversion via FFT. The application of this method in a
related model can be found in Chen (2007). The de Athayde and Flores (2002)
paper is the basis for the geometric properties of the higher moment tensors in
finance.\cr
}
\examples{
\dontrun{
data(dji30ret)
spec = gogarchspec(mean.model = list(demean = "constant"),
variance.model = list(model = "sGARCH", garchOrder = c(1,1), submodel = NULL),
distribution.model = list(distribution = "manig"), ica = "fastica")
fit = gogarchfit(spec = spec, data  = dji30ret[,1:4, drop = FALSE],
out.sample = 50, gfun = "tanh")
# The likelihood of the model
likelihood(fit)
# the GARCH coefficients of the independent factors
coef(fit)
# a news-impact surface plot
#ni = nisurface(fit, type = "cov", pair = c(1, 2), factor = c(1,2), plot = TRUE)
# the time varying correlation array
mc = rcor(fit)
# plot(mc[1,2,], type = "l")
# The moments of an equally weighted portfolio (subtract the out.sample from dimension)
gm = gportmoments(fit, weights = matrix(1/4, ncol = 4,
nrow = dim(dji30ret)[1]-50), debug = TRUE)
}}
\keyword{classes}
