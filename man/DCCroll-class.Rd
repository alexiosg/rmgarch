\name{DCCroll-class}
\docType{class}
\alias{DCCroll-class}
\alias{coef,DCCroll-method}
\alias{fitted,DCCroll-method}
\alias{likelihood,DCCroll-method}
\alias{plot,DCCroll,missing-method}
\alias{rcor,DCCroll-method}
\alias{rcov,DCCroll-method}
\alias{rshape,DCCroll-method}
\alias{rskew,DCCroll-method}
\alias{show,DCCroll-method}
\alias{sigma,DCCroll-method}

\title{class: DCC Roll Class}
\description{
The class is returned by calling the function \code{\link{dccroll}}.
}
\section{Slots}{
  \describe{
    \item{\code{mforecast}:}{Object of class \code{"vector"} Multivariate
    forecast list.}
    \item{\code{model}:}{Object of class \code{"vector"} Model specification
    list.}
  }
}
\section{Extends}{
Class \code{mGARCHroll}, directly.\\
Class \code{GARCHroll} object from the \pkg{rugarch} package, by class \code{mGARCHroll}, distance 2.\\
Class \code{rGARCH} object from the \pkg{rugarch} package, by class \code{mGARCHroll}, distance 3.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "DCCroll")}:
    The coefficient array across the rolling estimations with a T+0
    3rd dimension index label.}
    \item{fitted}{\code{signature(object = "DCCroll")}:
    The conditional mean forecast xts object (with the actual T+i forecast
    dates as index).}
    \item{likelihood}{\code{signature(object = "DCCroll")}:
    The log-likelihood across rolling estimations.}
    \item{plot}{\code{signature(x = "DCCroll", y = "missing")}:
    Plot method, given additional arguments \sQuote{series} and \sQuote{which}.}
    \item{rcor}{\code{signature(object = "DCCroll")}:
    The forecast dynamic conditional correlation array, with the T+i forecast
    index in the third dimension label. Optional argument \sQuote{type}
    determines whether to return \dQuote{R} for the correlation else will return
    the DCC Q matrix. A further argument \sQuote{output} allows to switch
    between \dQuote{array} and \dQuote{matrix} returned object.}
    \item{rcov}{\code{signature(object = "DCCroll")}:
    The forecast dynamic conditional covariance array, with the T+i forecast
    index in the 3rd dimension label. A further argument
    \sQuote{output} allows to switch between \dQuote{array}
    and \dQuote{matrix} returned object.}
  	\item{rshape}{\code{signature(object = "DCCroll")}:
    The multivariate distribution shape parameter(s).}
    \item{rskew}{\code{signature(object = "DCCroll")}:
    The multivariate distribution skew parameter(s). }
    \item{show}{\code{signature(object = "DCCroll")}:
    Summary. }
    \item{sigma}{\code{signature(object = "DCCroll")}:
    The conditional sigma forecast xts object (with the actual T+i forecast
    dates as index).}
	}
}
\references{
Engle, R.F. and Sheppard, K. 2001, Theoretical and empirical properties of
dynamic conditional correlation multivariate GARCH, \emph{NBER Working Paper}.\cr
}
\author{Alexios Galanos}
\keyword{classes}
