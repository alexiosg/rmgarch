\name{DCCforecast-class}
\docType{class}
\alias{DCCforecast-class}
\alias{fitted,DCCforecast-method}
\alias{sigma,DCCforecast-method}
\alias{rcor,DCCforecast-method}
\alias{rcov,DCCforecast-method}
\alias{rshape,DCCforecast-method}
\alias{rskew,DCCforecast-method}
\alias{show,DCCforecast-method}
\alias{plot,DCCforecast,missing-method}
\title{class: DCC Forecast Class}
\description{
The class is returned by calling the function \code{\link{dccforecast}}.
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
Class \code{"\linkS4class{mGARCHforecast}"}, directly.
Class \code{"\linkS4class{GARCHforecast}"}, by class "mGARCHforecast", distance 2.
Class \code{"\linkS4class{rGARCH}"}, by class "mGARCHforecast", distance 3.
}

\section{Methods}{
  \describe{
  	\item{rshape}{\code{signature(object = "DCCforecast")}:
    The multivariate distribution shape parameter(s).}
    \item{rskew}{\code{signature(object = "DCCforecast")}:
    The multivariate distribution skew parameter(s). }
    \item{fitted}{\code{signature(object = "DCCforecast")}:
    The conditional mean forecast array of dimensions n.ahead x n.assets
    by (n.roll+1). The thirds dimension of the array has the T+0 index label.}
    \item{sigma}{\code{signature(object = "DCCforecast")}:
    The conditional sigma forecast array of dimensions n.ahead x n.assets
    by (n.roll+1). The thirds dimension of the array has the T+0 index label.}
    \item{plot}{\code{signature(x = "DCCforecast", y = "missing")}:
    Plot method, given additional arguments \sQuote{series} and \sQuote{which}.}
    \item{rcor}{\code{signature(object = "DCCforecast")}:
    The forecast dynamic conditional correlation list of arrays of length
    (n.roll+1), with each array of dimensions n.assets x n.assets x n.ahead.
    The method takes on one additional argument \sQuote{type} (either \dQuote{R}
    for the correlation else will return the DCC Q matrix). A further argument
    \sQuote{output} allows to switch between \dQuote{array}
    and \dQuote{matrix} returned object.
    }
    \item{rcov}{\code{signature(object = "DCCforecast")}:
     The forecast dynamic conditional correlation list of arrays of length
    (n.roll+1), with each array of dimensions n.assets x n.assets x n.ahead.
    A further argument \sQuote{output} allows to switch between \dQuote{array}
    and \dQuote{matrix} returned object.
    }
    \item{show}{\code{signature(object = "DCCforecast")}:
    Summary. }
	 }
}
\references{
Engle, R.F. and Sheppard, K. 2001, Theoretical and empirical properties of
dynamic conditional correlation multivariate GARCH, \emph{NBER Working Paper}.\cr
}
\author{Alexios Galanos}
\keyword{classes}
