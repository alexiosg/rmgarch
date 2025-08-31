\name{fMoments-class}
\docType{class}
\alias{fMoments-class}
\alias{show,fMoments-method}
\alias{fitted,fMoments-method}
\alias{rcov,fMoments-method}
\alias{rcoskew,fMoments-method}
\alias{rcokurt,fMoments-method}
\title{Class \code{"fMoments"}}
\description{
Object returned from calling \code{\link{fmoments}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("fMoments", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{moments}:}{Object of class \code{"vector"} A list with the
    (roll+1) n-ahead forecast moment matrices.}
    \item{\code{model}:}{Object of class \code{"vector"} A list with details of the
    data-generating process.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "fMoments")}: Summary method. }
    \item{fitted}{\code{signature(object = "fMoments")}: Conditional mean
    forecast matrix. }
    \item{rcov}{\code{signature(object = "fMoments")}: Conditional covariance
    forecast array.}
    \item{rcoskew}{\code{signature(object = "fMoments")}: Conditional third
    co-moment array. }
    \item{rcokurt}{\code{signature(object = "fMoments")}: Conditional fourth
    co-moment array. }
    }
}
\author{Alexios Galanos}
\examples{
showClass("fMoments")
}
\keyword{classes}
