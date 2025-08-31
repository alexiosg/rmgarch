\name{goGARCHroll-class}
\docType{class}
\alias{goGARCHroll-class}
\alias{fitted,goGARCHroll-method}
\alias{sigma,goGARCHroll-method}
\alias{rcor,goGARCHroll-method}
\alias{rcov,goGARCHroll-method}
\alias{rcoskew,goGARCHroll-method}
\alias{rcokurt,goGARCHroll-method}
\alias{convolution,goGARCHroll-method}
\alias{gportmoments,goGARCHroll-method}
\alias{coef,goGARCHroll-method}
\title{class: GO-GARCH Roll Class}
\description{Class for the GO-GARCH Roll.}
\section{Objects from the Class}{
The class is returned by calling the function \code{\link{gogarchroll}}.
}
\section{Slots}{
  \describe{
    \item{\code{forecast}:}{Object of class \code{"vector"} which contains the
    rolling forecasts of the distributional parameters for each factor.}
    \item{\code{model}:}{Object of class \code{"vector"} containing details of the
    GOGARCH model specification.}
  }
}
\section{Extends}{
Class \code{mGARCHroll}, directly.\\
Class \code{GARCHroll} object from the \pkg{rugarch} package, by class \code{mGARCHroll}, distance 2.\\
Class \code{rGARCH} object from the \pkg{rugarch} package, by class \code{mGARCHroll}, distance 3.
}
\section{Methods}{
\describe{
    \item{coef}{\code{signature(object = "goGARCHroll")}:
    Extraction of independent factor GARCH model coefficients saved from the
    goGARCHfit objects (returns a list).}
  	\item{fitted}{\code{signature(object = "goGARCHroll")}:
  	Extracts the conditional fitted forecast values (returns an xts object with
  	index the actual forecast T+1 times).}
  	\item{sigma}{\code{signature(object = "goGARCHroll")}:
  	Extracts the conditional sigma forecast values (returns an xts object with
  	index the actual forecast T+1 times). Takes optional argument \dQuote{factors}
  	(default TRUE) denoting whether to return the factor conditional sigma or
  	the transformed sigma for the assets.}
    \item{rcov}{\code{signature(object = "goGARCHroll")}:
    Returns the time-varying n.asset x n.asset x (n.roll+1) covariance matrix in
    array format, where the third dimension labels are now the actual rolling
    n.ahead=1 forecast time indices (T+1). A further argument \sQuote{output} allows
    to switch between \dQuote{array} and \dQuote{matrix} returned object.}
    \item{rcor}{\code{signature(object = "goGARCHroll")}:
    Returns the time-varying n.asset x n.asset x (n.roll+1) correlation matrix
    in array format, where the third dimension labels are now the actual rolling
    n.ahead=1 forecast time indices (T+1). A further argument \sQuote{output}
    allows to switch between \dQuote{array} and \dQuote{matrix} returned object.}
    \item{rcoskew}{\code{signature(object = "goGARCHroll")}:
    Returns the time-varying n.asset x n.asset^2 x (n.roll+1) coskewness matrix
    in array format, where the third dimension labels are now the actual rolling
    n.ahead=1 forecast time indices (T+1). There is a \dQuote{standardize}
    option which indicates whether the coskewness should be standardized by the
    conditional sigma (see equations in vignette).}
    \item{rcokurt}{\code{signature(object = "goGARCHroll")}:
    Returns the time-varying n.asset x n.asset^3 x (n.roll+1) cokurtosis matrix
    in array format, where the third dimension labels are now the actual rolling
    n.ahead=1 forecast time indices (T+1). There is a \dQuote{standardize}
    option which indicates whether the cokurtosis should be standardized by the
    conditional sigma (see equations in vignette).}
    \item{gportmoments}{\code{signature(object = "goGARCHroll")}:\cr
    function:\cr
    \bold{gportmoments(object, weights)}\cr
    Calculates the first 4 standardized portfolio moments using the geometric
    properties of the model, given a matrix of asset weights with row dimension
    equal to the total rolling forecast horizon. Returns an xts object of
    dimensions (total rolling forecast) x 4 (moments), with the index denoting
    the T+1 actual forecast time. If the number of assets > 100, then the
    kurtosis is not returned (see cokurtosis restrictions below).}
    \item{convolution}{\code{signature(object = "goGARCHroll")}:\cr
    function:\cr
    \bold{convolution(object, weights, fft.step = 0.001, fft.by = 0.0001,
    fft.support = c(-1, 1), support.method = c("user", "adaptive"),
    use.ff = TRUE, cluster = NULL, trace = 0,...)}\cr
    The convolution method takes a goGARCHroll object and a weights vector or
    matrix and calculates the weighted density. If a vector is given, it must be
    the same length as the number of assets, otherwise a matrix with
    row dimension equal to the row dimension of total forecast horizon.
    In the case of the multivariate normal distribution, this simply returns the
    linear and quadratic transformation of the mean and covariance matrix, while
    in the multivariate affine NIG distribution this is based on the numerical
    inversion by FFT of the characteristic function. In that case, the
    \dQuote{fft.step} option determines the stepsize for tuning the
    characteristic function inversion, \dQuote{fft.by} determines the resolution
    for the equally spaced support given by \dQuote{fft.support}, while the use
    of the \dQuote{ff} package is recommended to avoid memory problems on some
    systems and is turned on via the \dQuote{use.ff} option.
    The \dQuote{support.method} option allows either a fixed support range to be
    given (option \sQuote{user}), else an adaptive method is used based on the
    min and max of the assets at each point in time at the 0.00001 and 1-0.00001
    quantiles. The range is equally spaced subject to the \dQuote{fft.by} value
    but the returned object no longer makes use of the \dQuote{ff} package returning
    instead a list. The option for parallel computation is available
    via the use of a cluster object as elsewhere in this package. Passing this
    object to the distribution methods (e.g. qfft) follows the same rules as
    the goGARCHforecast object, namely that the index is zero based.}
    \item{show}{\code{signature(object = "goGARCHroll")}: Summary. }
   }
}
\author{Alexios Galanos}
\keyword{classes}
