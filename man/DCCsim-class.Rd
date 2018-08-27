\name{DCCsim-class}
\docType{class}
\alias{DCCsim-class}
\alias{fitted,DCCsim-method}
\alias{rcor,DCCsim-method}
\alias{rcov,DCCsim-method}
\alias{sigma,DCCsim-method}
\alias{show,DCCsim-method}
\title{class: DCC Forecast Class}
\description{
The class is returned by calling the function \code{\link{dccsim}}.
}
\section{Slots}{
  \describe{
    \item{\code{msim}:}{Object of class \code{"vector"} Multivariate
    simulation list.}
    \item{\code{model}:}{Object of class \code{"vector"} Model specification
    list.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{mGARCHsim}"}, directly.
Class \code{"\linkS4class{GARCHsim}"}, by class "mGARCHsim", distance 2.
Class \code{"\linkS4class{rGARCH}"}, by class "mGARCHsim", distance 3.
}

\section{Methods}{
  \describe{
    \item{fitted}{\code{signature(object = "DCCsim")}:
    The conditional mean simulated data matrix given additional argument
    \sQuote{sim} denoting the simulation run (\code{m.sim}) to return values
    for.}
    \item{rcor}{\code{signature(object 	= "DCCsim")}:
    The simulated dynamic conditional correlation array given additional
    arguments \sQuote{sim} denoting the simulation run (\code{m.sim}) to
    return values for, and \sQuote{type} (either \dQuote{R} for the correlation
    else will return the Q matrix). A further argument \sQuote{output} allows to
    switch between \dQuote{array} and \dQuote{matrix} returned object.}
    \item{rcov}{\code{signature(object 	= "DCCsim")}:
    The simulated dynamic conditional covariance array given additional argument
    \sQuote{sim} denoting the simulation run (\code{m.sim}) to return values
    for. A further argument \sQuote{output} allows to switch between \dQuote{array}
    and \dQuote{matrix} returned object.}
    \item{sigma}{\code{signature(object = "DCCsim")}:
    The univariate simulated conditional sigma matrix given additional argument
    \sQuote{sim} (\code{m.sim}) denoting the simulation run to return values
    for.}
    \item{show}{\code{signature(object = "DCCsim")}:
    Summary. }
	}
}
\references{
Engle, R.F. and Sheppard, K. 2001, Theoretical and empirical properties of
dynamic conditional correlation multivariate GARCH, \emph{NBER Working Paper}.\cr
}
\author{Alexios Ghalanos}
\keyword{classes}
