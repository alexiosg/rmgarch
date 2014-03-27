\name{cGARCHsim-class}
\docType{class}
\alias{cGARCHsim-class}
\alias{fitted,cGARCHsim-method}
\alias{sigma,cGARCHsim-method}
\alias{rcor,cGARCHsim-method}
\alias{rcov,cGARCHsim-method}
\alias{show,cGARCHsim-method}
\title{class: Copula Simulation Class}
\description{
The class is returned by calling the function \code{\link{cgarchsim}}.
}
\section{Slots}{
  \describe{
	\item{\code{msim}:}{Object of class \code{"vector"} Multivariate simulation list.}
	\item{\code{model}:}{Object of class \code{"vector"} Model specification list. }
  }
}
\section{Extends}{
Class \code{"\linkS4class{mGARCHsim}"}, directly.
Class \code{"\linkS4class{GARCHsim}"}, by class "mGARCHsim", distance 2.
Class \code{"\linkS4class{rGARCH}"}, by class "mGARCHsim", distance 3.
}
\section{Methods}{
  \describe{
    \item{fitted}{\code{signature(object = "cGARCHsim")}: 
    The simulated conditional returns matrix given. Takes optional argument 
    \dQuote{sim} indicating the simulation run to return (from the m.sim option 
    of the \code{\link{cgarchsim}} method.}
    \item{sigma}{\code{signature(object = "cGARCHfit")}: 
    The simulated conditional sigma matrix given. Takes optional argument 
    \dQuote{sim} indicating the simulation run to return (from the m.sim option 
    of the \code{\link{cgarchsim}} method.}
    \item{rcor}{\code{signature(object = "cGARCHsim")}: 
    The simulated conditional correlation array (for DCC type). Takes optional 
    argument \dQuote{sim} indicating the simulation run to return (from the 
    m.sim option of the \code{\link{cgarchsim}} method}
    \item{rcov}{\code{signature(object = "cGARCHsim")}: 
    The simulated conditional covariance array. Takes optional argument 
    \dQuote{sim} indicating the simulation run to return (from the m.sim option 
    of the \code{\link{cgarchsim}} method.}
    \item{show}{\code{signature(object = "cGARCHsim")}: 
    Summary. }
	}
}
\references{
Joe, H. \emph{Multivariate Models and Dependence Concepts}, 1997, 
Chapman \& Hall, London.\cr
Genest, C., Ghoudi, K. and Rivest, L. \emph{A semiparametric estimation 
procedure of dependence parameters in multivariate families of distributions}, 
1995, Biometrika, 82, 543-552.\cr
}
\author{Alexios Ghalanos}
\keyword{classes}