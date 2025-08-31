\name{goGARCHforecast-class}
\docType{class}
\alias{goGARCHforecast-class}
\alias{coef,goGARCHforecast-method}
\alias{fitted,goGARCHforecast-method}
\alias{sigma,goGARCHforecast-method}
\alias{as.matrix,goGARCHforecast-method}
\alias{convolution,goGARCHforecast-method}
\alias{gportmoments,goGARCHforecast-method}
\alias{rcokurt,goGARCHforecast-method}
\alias{rcor,goGARCHforecast-method}
\alias{rcoskew,goGARCHforecast-method}
\alias{rcov,goGARCHforecast-method}
\alias{betacovar,goGARCHforecast-method}
\alias{betacoskew,goGARCHforecast-method}
\alias{betacokurt,goGARCHforecast-method}
\alias{show,goGARCHforecast-method}
\title{class: GO-GARCH Forecast Class}
\description{Class for the GO-GARCH forecast.}
\section{Objects from the Class}{
The class is returned by calling the function \code{\link{gogarchforecast}}.
}
\section{Slots}{
  \describe{
    \item{\code{mforecast}:}{Multivariate forecast object.}
    \item{\code{model}:}{Object of class \code{"vector"} containing details of the
    GOGARCH model specification.}
  }
}
\section{Extends}{
Class \code{mGARCHforecast}, directly.\\
Class \code{GARCHforecast} object from the \pkg{rugarch} package, by class \code{mGARCHforecast}, distance 2.\\
Class \code{rGARCH} object from the \pkg{rugarch} package, by class \code{mGARCHforecast}, distance 3.
}
\section{Methods}{
  \describe{
    \item{convolution}{\code{signature(object = "goGARCHforecast")}:\cr
    function:\cr
    \bold{convolution(object, weights, fft.step = 0.001, fft.by = 0.0001,
    fft.support = c(-1, 1), support.method = c("user", "adaptive"),
    use.ff = TRUE, cluster = NULL, trace = 0,...)}\cr
    The convolution method takes a goGARCHforecast object and a weights vector
    or matrix and calculates the weighted density. If a vector is given, it must
    be the same length as the number of assets, otherwise a matrix with
    row dimension equal to the total forecast horizon.
    In the case of the multivariate normal distribution, this simply returns the
    linear and quadratic transformation of the mean and covariance matrix,
    while in the multivariate affine NIG distribution this is based on the
    numerical inversion by FFT of the characteristic function. In that case,
    the \dQuote{fft.step} option determines the stepsize for tuning the
    characteristic function inversion, \dQuote{fft.by} determines the resolution
    for the equally spaced support given by \dQuote{fft.support}, while the use
    of the \dQuote{ff} package is recommended to avoid memory problems on some
    systems and is turned on via the \dQuote{use.ff} option. The \dQuote{support.method}
    option allows either a fixed support range to be given (option \sQuote{user}),
    else an adaptive method is used based on the min and max of the assets at
    each point in time at the 0.00001 and 1-0.00001 quantiles. The range is
    equally spaced subject to the \dQuote{fft.by} value but the returned object
    no longer makes of the \dQuote{ff} package returning instead a list.
    The option for parallel computation is available via the use of a
    cluster object as elsewhere in this package. There is no special treatment of
    the forecast type here (unconditional or rolling), since either n.ahead with
    no roll or rolling with 1-ahead only choices are available for the
    \code{gogarchforecast} method. This means that the stored object does not
    distringuish between an unconditional or rolling forecast, calculating the
    density for all points (see note).
    }
    \item{gportmoments}{\code{signature(object = "goGARCHforecast")}:\cr
    function:\cr
    \bold{gportmoments(object, weights)}\cr
    Calculates the first 4 standardized portfolio moments using the geometric
    properties of the model, given a matrix of asset weights with row dimension
    equal to the forecast n.ahead or n.roll horizon. Returns an array of
    dimensions n.ahead x 4 (moments) x n.roll, with the third array dimension
    labelled with the T+0 index times. If the number of assets > 100, then the
    kurtosis is not returned (see cokurtosis restrictions below).}
    \item{rcoskew}{\code{signature(object = "goGARCHforecast")}:\cr
    function:\cr
    \bold{rcoskew(object, standardize = TRUE, from = 1, to = 1, roll = 0)}\cr
    Returns the 'time-varying'  NxN^2 (coskewness tensor) x (to:from|roll) in
    array format. The \dQuote{from} and \dQuote{to} options indicate the time
    indices for which to return the array and \dQuote{roll} the rolling index
    (base=0). The third dimension array label denotes the T+i (i=from:to)
    forecast horizon given the T+0 roll index which is returned as an attribute
    (attr(,\dQuote{T+0})) of the array. The \dQuote{standardize} option
    indicates whether the coskewness should be standardized by the
    conditional sigma (see equations in vignette). It is also possible to set
    roll to the character \sQuote{all} in which case all the rolling 1-ahead
    forecasts are returned in an n by n^2 by (n.roll+1) array with 3rd
    dimension label the T+0 dates (instead of being an attribute).}
    \item{rcokurt}{\code{signature(object = "goGARCHforecast")}:\cr
    function:\cr
    \bold{rcokurt(object, standardize = TRUE, from = 1, to = 1, roll = 0)}\cr
    Returns the 'time-varying'  NxN^3 (cokurtosis tensor) x (to:from|roll) in
    array format. The \dQuote{from} and \dQuote{to} options indicate the time
    indices for which to return the array and \dQuote{roll} the rolling index
    (base=0).  Because of memory issues, this is only returned when the number
    of assets are less than 100. The third dimension array label denotes the
    T+i (i=from:to) forecast horizon given the T+0 roll index which is returned
    as an attribute (attr(,\dQuote{T+0})) of the array. The \dQuote{standardize}
    option indicates whether the cokurtosis should be standardized by the
    conditional sigma (see equations in vignette). It is also possible to set
    roll to the character \sQuote{all} in which case all the rolling 1-ahead
    forecasts are returned in an n by n^3 by (n.roll+1) array with 3rd
    dimension label the T+0 dates (instead of being an attribute.}
    \item{rcov}{\code{signature(object = "goGARCHforecast")}:
    Returns the conditional covariances, in a list of length (n.roll+1), with
    names the T+0 index, and each list slot having an array of dimensions
    n.asset x n.asset x n.ahead, with the third array dimension labelled
    as T+i (i>0). A further argument \sQuote{output} allows to switch
    between \dQuote{array} and \dQuote{matrix} returned object.}
    \item{rcor}{\code{signature(object = "goGARCHforecast")}:
    Returns the conditional correlations, in a list of length (n.roll+1), with
    names the T+0 index, and each list slot having an array of dimensions
    n.asset x n.asset x n.ahead, with the third array dimension labelled
    as T+i (i>0). A further argument \sQuote{output} allows to switch
    between \dQuote{array} and \dQuote{matrix} returned object.}
    \item{coef}{\code{signature(object = "goGARCHforecast")}:
	Extraction of independent factor GARCH model coefficients saved from the
	goGARCHfit object.}
	\item{fitted}{\code{signature(object = "goGARCHforecast")}:
	Extracts the conditional mean forecast values. Returns an n.ahead
	x n.assets x (n.roll+1) array where the third dimension array labels are the
	T+0 index times.}
	\item{sigma}{\code{signature(object = "goGARCHforecast")}:
	Extracts the conditional sigma forecast values. Returns an n.ahead
	x n.assets x (n.roll+1) array where the third dimension array labels are the
	T+0 index times. Takes optional argument \dQuote{factors} (default TRUE)
	denoting whether to return the factor conditional sigma or the transformed
	sigma for the assets.}
	\item{as.matrix}{\code{signature(x = "goGARCHforecast")}:\cr
	function:\cr
	\bold{as.matrix(x, which = "A")}\cr
	This returns four types of matrices relating to the estimation of the
	independent components in the GO-GARCH model. Valid choices are \dQuote{A}
	for the mixing matrix, \dQuote{W} for the unmixing matrix, \dQuote{U} for the
	rotational matrix and \dQuote{K} for the whitening matrix, \dQuote{Kinv} for
	the de-whitening matrix.}
	\item{betacovar}{\code{signature(object = "goGARCHforecast")}:
	function:\cr
	\bold{betacovar(object, weights, asset = 1)}\cr
	Returns the covariance beta given a matrix (of length equal to the number of
	rows of the forecast horizon, or vector which is then recycled to the number
	of rows of the forecast horizon) of benchmark weights and the asset number.}
	\item{betacoskew}{\code{signature(object = "goGARCHforecast")}:
	function:\cr
	\bold{betacoskew(object, weights, asset = 1)}\cr
	Returns the coskewness beta given a matrix (of length equal to the number of
	rows of the forecast horizon, or vector which is then recycled to the number
	of rows of the forecast horizon) of benchmark weights and the asset number.}
	\item{betacokurt}{\code{signature(object = "goGARCHforecast")}:
	function:\cr
	\bold{betacokurt(object, weights, asset = 1)}\cr
	Returns the cokurtosis beta given a matrix (of length equal to the number of
	rows of the forecast horizon, or vector which is then recycled to the number
	of rows of the forecast horizon) of benchmark weights and the asset number.}
    \item{show}{\code{signature(object = "goGARCHforecast")}: Summary method. }
	}
}
\references{
Chen, Y., Hardle, W., and Spokoiny, V. 2010, GHICA-Risk analysis with GH
distributions and independent components, \emph{Journal of Empirical Finance},
\bold{17(2)}, 255--269.\cr
de Athayde, G.M. and Flores Jr, R.G. 2002, On Certain Geometric Aspects of
Portfolio Optimisation with Higher Moments, \emph{mimeo}.\cr
Ghalanos, A., Rossi, E., and Urga, G. 2015, Independent factor autoregressive
conditional density model, \emph{Econometric Reviews}, \bold{34}(5), 594--616.\cr
Paolella, M.S. 2007, Intermediate Probability - A Computational Approach,
\emph{Wiley-Interscience}.\cr
}
\author{Alexios Galanos}
\note{
The reference by Chen et al (2010) and Paolella (2007) contains more details on
the algorithm for the characteristic function inversion via FFT. The de Athayde
and Flores (2002) paper is the basis for some of the geometric properties of
the higher moment tensors. The paper by Ghalanos et al (2013) contains more
specific details.\cr
Forecasts are carried out on the time varying parameters of the factor
distributions, and then scaled and transformed to those of the assets after
adding back the mean forecast (which is either a constant or the AR/VAR mean
forecast).
}
\keyword{classes}
