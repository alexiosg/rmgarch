\name{goGARCHsim-class}
\docType{class}
\alias{goGARCHsim-class}
\alias{convolution,goGARCHsim-method}
\alias{gportmoments,goGARCHsim-method}
\alias{rcor,goGARCHsim-method}
\alias{rcoskew,goGARCHsim-method}
\alias{rcokurt,goGARCHsim-method}
\alias{rcov,goGARCHsim-method}
\alias{as.matrix,goGARCHsim-method}
\title{class: GO-GARCH Simulation Class}
\description{Class for the GO-GARCH Simulation.}
\section{Objects from the Class}{
The class is returned by calling the function \code{\link{gogarchsim}}.
}
\section{Slots}{
  \describe{
    \item{\code{msim}:}{Object of class \code{"vector"} The multivariate
    simulation list.}
    \item{\code{model}:}{Object of class \code{"vector"} containing details of the
    GOGARCH model specification.}
  }
}
\section{Extends}{
Class \code{mGARCHsim}, directly.\\
Class \code{GARCHsim} object from the \pkg{rugarch} package, by class \code{mGARCHsim}, distance 2.\\
Class \code{rGARCH} object from the \pkg{rugarch} package, by class \code{mGARCHsim}, distance 3.
}
\section{Methods}{
\describe{
\item{convolution}{\code{signature(object = "goGARCHsim")}:\cr
	function:\cr
    \bold{convolution(object, weights, fft.step = 0.001, fft.by = 0.0001,
    fft.support = c(-1, 1), support.method = c("user", "adaptive"),
    use.ff = TRUE, sim = 1, cluster = NULL, trace = 0,...)}\cr
    The convolution method takes a goGARCHsim object and a weights vector and
    calculates the weighted density. The vector must be the same length as the
    number of assets.
    The \dQuote{sim} option indicates the simulation index to  use, given
    the \dQuote{m.sim} option chosen in the call to the simulation function.
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
    equally spaced subject to the \dQuote{fft.by} value but the returned object no
    longer makes use of the \dQuote{ff} package returning instead a list.
    Finally, the option for parallel computation is available via the use of a
    cluster object as elsewhere in this package.}
\item{gportmoments}{\code{signature(object = "goGARCHsim")}:\cr
    function:\cr
    \bold{gportmoments(object, weights, sim = 1)}\cr
    Calculates the first 3 portfolio moments using the geometric properties of
    the model, given a matrix of asset weights with row dimension equal to the
    row dimension of the filtered dataset (i.e. less any lags).
    The \dQuote{sim} option indicates the simulation index to use,
    given the \dQuote{m.sim} option chosen in the call to the simulation function.
    }
\item{rcoskew}{\code{signature(object = "goGARCHsim")}:\cr
    function:\cr
    \bold{rcoskew(object, from = 1, to = 1, sim = 1)}\cr
    Returns the 'time-varying'  NxN^2 coskewness tensor in array format.
    The \dQuote{from} and \dQuote{to} options indicate the time indices for
    which to return the arrays. Because of memory issues, this is limited to 100
    indices.  The \dQuote{sim} option indicates the simulation index to use,
    given the \dQuote{m.sim} option chosen in the call to the simulation
    function.}
\item{rcokurt}{\code{signature(object = "goGARCHsim")}:\cr
    function:\cr
    \bold{rcokurt(object, standardize = TRUE, from = 1, to = 1)}\cr
    Returns the 'time-varying'  NxN^3 cokurtosis tensor in array format. The
    \dQuote{from} and \dQuote{to} options indicate the time indices for which
    to return the arrays. Because of memory issues, this is limited to models
	with less than 20 assets.}
\item{rcov}{\code{signature(object = "goGARCHsim")}:
	Returns the time-varying NxN covariance matrix in array format.  There is an
	additional \dQuote{sim} option which indicates the simulation index to
    use, given the \dQuote{m.sim} option chosen in the call to the simulation
    function. A further argument \sQuote{output} allows to switch between
    \dQuote{array} and \dQuote{matrix} returned object.}
\item{rcor}{\code{signature(object = "goGARCHsim")}:
	Returns the time-varying NxN correlation matrix in array format.
	There is an additional \dQuote{sim} option which indicates the simulation
	index to use, given the \dQuote{m.sim} option chosen in the call to the
	simulation function. A further argument \sQuote{output} allows to switch
	between \dQuote{array} and \dQuote{matrix} returned object.}
\item{as.matrix}{\code{signature(x = "goGARCHsim")}:\cr
	function:\cr
	\bold{as.matrix(x, which = "A")}\cr
	This returns four types of matrices relating to the estimation of the
	independent components in the GO-GARCH model. Valid choices are \dQuote{A}
	for the mixing matrix, \dQuote{W} for the unmixing matrix, \dQuote{U} for the
	rotational matrix and \dQuote{K} for the whitening matrix, \dQuote{Kinv} for
	the de-whitening matrix.}
	}
}
\author{Alexios Galanos}
\keyword{classes}
