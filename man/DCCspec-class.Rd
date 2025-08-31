\name{DCCspec-class}
\docType{class}
\alias{DCCspec-class}
\alias{show,DCCspec-method}
\alias{setfixed<-,DCCspec,vector-method}
\alias{setstart<-,DCCspec,vector-method}
\title{class: DCC Specification Class}
\description{
The class is returned by calling the function \code{\link{dccspec}}.
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{Object of class \code{"vector"} The multivariate model
    specification list. }
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
    \item{setfixed<-}{\code{signature(object = "DCCspec", value = "vector")}:
    Set fixed second stage parameters.}
    \item{setstart<-}{\code{signature(object = "DCCspec", value = "vector")}:
    Set starting second stage parameters.}
    \item{show}{\code{signature(object = "DCCspec")}:
    Summary.}
	 }
}
\references{
Croux, C. and Joossens, K. 2008, Robust estimation of the vector autoregressive
model by a least trimmed squares procedure, \emph{COMPSTAT}, 489--501.\cr
Cappiello, L., Engle, R.F. and Sheppard, K. 2006, Asymmetric dynamics in the
correlations of global equity and bond returns, \emph{Journal of Financial
Econometrics} \bold{4}, 537--572.\cr
Engle, R.F. and Sheppard, K. 2001, Theoretical and empirical properties of
dynamic conditional correlation multivariate GARCH, \emph{NBER Working Paper}.\cr
}
\note{
The \sQuote{umodel} list is absorbed into the \sQuote{model} list in all other
DCC classes.
}
\author{Alexios Galanos}
\keyword{classes}
