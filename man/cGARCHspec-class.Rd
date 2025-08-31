\name{cGARCHspec-class}
\docType{class}
\alias{cGARCHspec-class}
\alias{show,cGARCHspec-method}
\alias{setfixed<-,cGARCHspec,vector-method}
\alias{setstart<-,cGARCHspec,vector-method}
\title{class: Copula Specification Class}
\description{
The class is returned by calling the function \code{\link{cgarchspec}}.
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"vector"}
    The multivariate model specification.}
    \item{\code{umodel}:}{Object of class \code{\link[rugarch:uGARCHmultispec-class]{uGARCHmultispec}}
    The univariate model specification.}
  }
}
\section{Extends}{
Class \code{\linkS4class{mGARCHspec}}, directly.
Class \code{\link[rugarch:GARCHspec-class]{GARCHspec}}, by class "mGARCHspec", distance 2.
Class \code{\link[rugarch:rGARCH-class]{rGARCH}}, by class "mGARCHspec", distance 3.
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "cGARCHspec")}:
    Summary.}
    \item{setfixed<-}{\code{signature(object = "cGARCHspec", value = "vector")}:
    Set fixed second stage parameters.}
    \item{setstart<-}{\code{signature(object = "cGARCHspec", value = "vector")}:
    Set starting second stage parameters.}
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
\keyword{classes}
