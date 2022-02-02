\name{cGARCHfilter-class}
\docType{class}
\alias{cGARCHfilter-class}
\alias{coef,cGARCHfilter-method}
\alias{fitted,cGARCHfilter-method}
\alias{likelihood,cGARCHfilter-method}
\alias{rcor,cGARCHfilter-method}
\alias{rcov,cGARCHfilter-method}
\alias{rskew,cGARCHfilter-method}
\alias{rshape,cGARCHfilter-method}
\alias{residuals,cGARCHfilter-method}
\alias{show,cGARCHfilter-method}
\alias{sigma,cGARCHfilter-method}
\title{class: Copula Filter Class}
\description{
The class is returned by calling the function \code{\link{cgarchfilter}}.
}
\section{Slots}{
  \describe{
    \item{\code{mfilter}:}{Object of class \code{"vector"} Multivariate 
    filter list. }
    \item{\code{model}:}{Object of class \code{"vector"} Model specification 
    list. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{mGARCHfilter}"}, directly.
Class \code{"\linkS4class{GARCHfilter}"}, by class "mGARCHfilter", distance 2.
Class \code{"\linkS4class{rGARCH}"}, by class "mGARCHfilter", distance 3.
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "cGARCHfilter")}: 
    The coefficient vector (see note). }
    \item{fitted}{\code{signature(object = "cGARCHfilter")}:  
    The conditional mean filtered data (xts object). }
    \item{likelihood}{\code{signature(object = "cGARCHfilter")}: 
    The joint likelihood. }
    \item{rcor}{\code{signature(object = "cGARCHfilter")}: 
    The conditional correlation array with third dimension labels the time index.}
    \item{rcov}{\code{signature(object = "cGARCHfilter")}: 
    The conditional covariance array with third dimension labels the time index.}
    \item{residuals}{\code{signature(object = "cGARCHfilter")}: 
    The model residuals (xts object).}
    \item{show}{\code{signature(object = "cGARCHfilter")}: 
    Summary. }
    \item{sigma}{\code{signature(object = "cGARCHfilter")}: 
    The model conditional sigma (xts object). }
    \item{rshape}{\code{signature(object = "cGARCHfilter")}: 
    The multivariate distribution shape parameter(s).}
    \item{rskew}{\code{signature(object = "cGARCHfilter")}: 
    The multivariate distribution skew parameter(s). }
	}
}
\references{
Joe, H. \emph{Multivariate Models and Dependence Concepts}, 1997, 
Chapman \& Hall, London.\cr
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
