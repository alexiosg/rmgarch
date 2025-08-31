\name{DCCfilter-class}
\docType{class}
\alias{DCCfilter-class}
\alias{coef,DCCfilter-method}
\alias{fitted,DCCfilter-method}
\alias{sigma,DCCfilter-method}
\alias{likelihood,DCCfilter-method}
\alias{nisurface,DCCfilter-method}
\alias{plot,DCCfilter,missing-method}
\alias{rcor,DCCfilter-method}
\alias{rcov,DCCfilter-method}
\alias{residuals,DCCfilter-method}
\alias{rshape,DCCfilter-method}
\alias{rskew,DCCfilter-method}
\alias{show,DCCfilter-method}
\title{class: DCC Filter Class}
\description{
The class is returned by calling the function \code{\link{dccfilter}}.
}
\section{Slots}{
  \describe{
    \item{\code{mfilter}:}{Object of class \code{"vector"}. Multivariate filter
    list.}
    \item{\code{model}:}{Object of class \code{"vector"}. Model specification
    list.}
  }
}
\section{Extends}{
Class \code{\link[rugarch:S4classes]{mGARCHfilter}}, directly.\\
Class \code{\link[rugarch:S4classes]{GARCHfilter}}, by class "mGARCHfilter", distance 2.\\
Class \code{\link[rugarch:S4classes]{rGARCH}}, by class "mGARCHfilter", distance 3.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "DCCfilter")}
    The coefficient vector (see note). }
    \item{likelihood}{\code{signature(object = "DCCfilter")}:
    The joint likelihood. }
    \item{rshape}{\code{signature(object = "DCCfilter")}:
    The multivariate distribution shape parameter(s).}
    \item{rskew}{\code{signature(object = "DCCfilter")}:
    The multivariate distribution skew parameter(s). }
    \item{fitted}{\code{signature(object = "DCCfilter")}:
    The filtered conditional mean xts object. }
    \item{sigma}{\code{signature(object = "DCCfilter")}:
    The filtered conditional sigma xts object. }
    \item{residuals}{\code{signature(object = "DCCfilter")}:
    The filtered conditional mean residuals xts object. }
    \item{plot}{\code{signature(x = "DCCfilter", y = "missing")}:
    Plot method, given additional arguments \sQuote{series} and \sQuote{which}.}
    \item{infocriteria}{\code{signature(object = "DCCfilter")}:
    Information criteria. }
    \item{rcor}{\code{signature(object = "DCCfilter")}:
    The filtered dynamic conditional correlation array given additional
    argument \sQuote{type} (either \dQuote{R} for the
    correlation else will return the \dQuote{Q} matrix). The third dimension
    label of the array gives the time index (from which it is then possible to
    construct pairwise xts objects for example). A further argument \sQuote{output}
    allows to switch between \dQuote{array} and \dQuote{matrix} returned object.}
    \item{rcov}{\code{signature(object = "DCCfilter")}:
    The filtered dynamic conditional covariance array. The third dimension
    label of the array gives the time index (from which it is then possible to
    construct pairwise xts objects for example). A further argument \sQuote{output}
    allows to switch between \dQuote{array} and \dQuote{matrix} returned object.}
    \item{show}{\code{signature(object = "DCCfilter")}:
    Summary. }
    \item{nisurface}{\code{signature(object = "DCCfilter")}:
    The news impact surface plot given additional arguments \sQuote{type} with
    either \dQuote{cov} or \dQuote{cor} (for the covariance and correlation news
    impact respectively), \sQuote{pair} denoting the asset pair (defaults to
    c(1,2)), \sQuote{plot} (logical) and \sQuote{plot.type} with a choice of
    either \dQuote{surface} or \dQuote{contour}. }
	 }
}
\references{
Engle, R.F. and Sheppard, K. 2001, Theoretical and empirical properties of
dynamic conditional correlation multivariate GARCH, \emph{NBER Working Paper}.\cr
}
\note{
The \sQuote{coef} method takes additional argument \sQuote{type} with valid
values \sQuote{garch} for the univariate garch parameters,  \sQuote{dcc} for the
second stage dcc parameters and by default returns all the parameters in a named
vector.
}
\author{Alexios Galanos}
\keyword{classes}
