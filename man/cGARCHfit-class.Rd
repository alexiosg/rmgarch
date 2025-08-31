\name{cGARCHfit-class}
\docType{class}
\alias{cGARCHfit-class}
\alias{coef,cGARCHfit-method}
\alias{fitted,cGARCHfit-method}
\alias{likelihood,cGARCHfit-method}
\alias{rcor,cGARCHfit-method}
\alias{rcov,cGARCHfit-method}
\alias{rskew,cGARCHfit-method}
\alias{rshape,cGARCHfit-method}
\alias{residuals,cGARCHfit-method}
\alias{show,cGARCHfit-method}
\alias{sigma,cGARCHfit-method}
\title{class: Copula Fit Class}
\description{
The class is returned by calling the function \code{\link{cgarchfit}}.
}
\section{Slots}{
  \describe{
    \item{\code{mfit}:}{Object of class \code{"vector"} Multivariate fit list. }
    \item{\code{model}:}{Object of class \code{"vector"} Model specification list. }
  }
}
\section{Extends}{
Class \code{mGARCHfit}, directly.\\
Class \code{GARCHfit} object from the \pkg{rugarch} package, by class \code{mGARCHfit}, distance 2.\\
Class \code{rGARCH} object from the \pkg{rugarch} package, by class \code{mGARCHfit}, distance 3.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "cGARCHfit")}:
    The coefficient vector (see note). }
    \item{fitted}{\code{signature(object = "cGARCHfit")}:
    The conditional mean fitted data (xts object). }
    \item{likelihood}{\code{signature(object = "cGARCHfit")}:
    The joint likelihood. }
    \item{rcor}{\code{signature(object = "cGARCHfit")}:
    The conditional correlation array with third dimension labels the time index.
    A further argument \sQuote{output} allows to switch between \dQuote{array}
    and \dQuote{matrix} returned object.}
    \item{rcov}{\code{signature(object = "cGARCHfit")}:
    The conditional covariance array with third dimension labels the time index.
    A further argument \sQuote{output} allows to switch between \dQuote{array}
    and \dQuote{matrix} returned object.}
    \item{rshape}{\code{signature(object = "cGARCHfit")}:
    The multivariate distribution shape parameter(s).}
    \item{rskew}{\code{signature(object = "cGARCHfit")}:
    The multivariate distribution skew parameter(s). }
    \item{residuals}{\code{signature(object = "cGARCHfit")}:
    The model residuals (xts object).}
    \item{show}{\code{signature(object = "cGARCHfit")}:
    Summary.}
    \item{sigma}{\code{signature(object = "cGARCHfit")}:
    The model conditional sigma (xts object). }
    }
}
\references{
Joe, H. \emph{Multivariate Models and Dependence Concepts}, 1997,
Chapman & Hall, London.\cr
Genest, C., Ghoudi, K. and Rivest, L. \emph{A semiparametric estimation
procedure of dependence parameters in multivariate families of distributions},
1995, Biometrika, 82, 543-552.\cr
}
\author{Alexios Galanos}
\note{
The \sQuote{coef} method takes additional argument \sQuote{type} with valid
values \sQuote{garch} for the garch parameters, \sQuote{dcc} for the second
stage parameters and by default returns all the parameters in a named vector.
}
\keyword{classes}
