\name{DCCfit-class}
\docType{class}
\alias{DCCfit-class}
\alias{coef,DCCfit-method}
\alias{fitted,DCCfit-method}
\alias{sigma,DCCfit-method}
\alias{infocriteria,DCCfit-method}
\alias{likelihood,DCCfit-method}
\alias{nisurface,DCCfit-method}
\alias{plot,DCCfit,missing-method}
\alias{rcor,DCCfit-method}
\alias{rcov,DCCfit-method}
\alias{residuals,DCCfit-method}
\alias{rshape,DCCfit-method}
\alias{rshape}
\alias{rskew,DCCfit-method}
\alias{rskew}
\alias{show,DCCfit-method}
\title{class: DCC Fit Class}
\description{
The class is returned by calling the function \code{\link{dccfit}}.
}
\section{Slots}{
  \describe{
    \item{\code{mfit}:}{Object of class \code{"vector"} Multivariate filter
    list.}
    \item{\code{model}:}{Object of class \code{"vector"}  Model specification
    list. }
  }
}
\section{Extends}{
Class \code{\linkS4class{mGARCHfit}}, directly.
Class \code{\link[rugarch:GARCHfit-class]{GARCHfit}}, by class "mGARCHfit", distance 2.
Class \code{\link[rugarch:rGARCH-class]{rGARCH}}, by class "mGARCHfit", distance 3.
}
\section{Methods}{
  \describe{
   \item{coef}{\code{signature(object = "DCCfit")}:
    The coefficient vector (see note). }
    \item{likelihood}{\code{signature(object = "DCCfit")}:
    The joint likelihood. }
    \item{rshape}{\code{signature(object = "DCCfit")}:
    The multivariate distribution shape parameter(s).}
    \item{rskew}{\code{signature(object = "DCCfit")}:
    The multivariate distribution skew parameter(s). }
    \item{fitted}{\code{signature(object = "DCCfit")}:
    The fitted conditional mean xts object. }
    \item{sigma}{\code{signature(object = "DCCfit")}:
    The fitted conditional GARCH sigma xts object. }
    \item{residuals}{\code{signature(object = "DCCfit")}:
    The fitted conditional mean residuals xts object. }
    \item{plot}{\code{signature(x = "DCCfit", y = "missing")}:
    Plot method, given additional arguments \sQuote{series} and \sQuote{which}.}
    \item{infocriteria}{\code{signature(object = "DCCfit")}:
    Information criteria. }
    \item{rcor}{\code{signature(object = "DCCfit")}:
    The fitted dynamic conditional correlation array given additional
    arguments \sQuote{type} (either \dQuote{R} for the
    correlation else will return the Q matrix). The third dimension label of the
    array gives the time index (from which it is then possible to construct
    pairwise xts objects for example). The argument \sQuote{output} can be either
    \dQuote{array} (default) or \dQuote{matrix} in which case the array is flattened
    and the lower diagonal time varying values are returned (and if a date exists,
    then the returned object is of class xts).
    }
    \item{rcov}{\code{signature(object = "DCCfit")}:
    The fitted dynamic conditional covariance array. The third dimension label
    of the array gives the time index (from which it is then possible to
    construct pairwise xts objects for example). The argument \sQuote{output} can be either
    \dQuote{array} (default) or \dQuote{matrix} in which case the array is flattened
    and the lower and main diagonal time varying values are returned (and if a date exists,
    then the returned object is of class xts).
    }
    \item{show}{\code{signature(object = "DCCfit")}:
    Summary. }
    \item{nisurface}{\code{signature(object = "DCCfit")}:
    The news impact surface plot given additional arguments \sQuote{type} with
    either \dQuote{cov} or \dQuote{cor} (for the covariance and correlation news
    impact respectively), \sQuote{pair} (defaults to c(1,2)), \sQuote{plot}
    (logical) and \sQuote{plot.type} with a choice of either \dQuote{surface} or
    \dQuote{contour}. }
    }
}
\references{
Engle, R.F. and Sheppard, K. 2001, Theoretical and empirical properties of
dynamic conditional correlation multivariate GARCH, \emph{NBER Working Paper}.\cr
}
\author{Alexios Galanos}
\note{
The \sQuote{coef} method takes additional argument \sQuote{type} with valid
values \sQuote{garch} for the univariate garch parameters,  \sQuote{dcc} for the
second stage dcc parameters and by default returns all the parameters in a named
vector.
}
\keyword{classes}
