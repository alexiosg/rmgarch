\name{goGARCHspec-class}
\docType{class}
\alias{goGARCHspec-class}
\alias{show,goGARCHspec-method}
\title{class: GO-GARCH Specification Class}
\description{Class for the GO-GARCH specification.}
\section{Objects from the Class}{
The class is returned by calling the function \code{\linkS4class{goGARCHspec}}.
}
\section{Slots}{
  \describe{
    \item{\code{model}:}{ Multivariate model specification.}
    \item{\code{umodel}:}{ Univariate model specification.}
  }
}
\section{Extends}{
Class \code{"\linkS4class{mGARCHspec}"}, directly.
Class \code{"\linkS4class{GARCHspec}"}, by class "mGARCHspec", distance 2.
Class \code{"\linkS4class{rGARCH}"}, by class "mGARCHspec", distance 3.
}
\section{Methods}{
\describe{
\item{show}{\code{signature(object = "goGARCHspec")}: Summary method. }
}
}
\references{
van der Weide, R. 2002, GO-GARCH: a multivariate generalized orthogonal GARCH 
model, \emph{Journal of Applied Econometrics},  549--564.\cr
Zhang, K. and Chan, L. 2009, Efficient factor GARCH models and factor-DCC models, 
\emph{Quantitative Finance}, 71--91.\cr
Broda, S.A. and Paolella, M.S. 2009, CHICAGO: A Fast and Accurate Method for 
Portfolio Risk Calculation, \emph{Journal of Financial Econometrics}, 412--436.\cr
Ghalanos, A. and Rossi, E. and Urga, G. 2011, Independent Factor Autoregressive 
Conditional Density Model, \emph{Pending--submitted}.
}
\author{Alexios Galanos}
\note{
The mixing matrix in the GO-GARCH model implemented in the rmgarch package is 
based on non-parametric independent component analysis (ICA) methodology. 
The estimation is a 2-stage methodology described in Broda and Paolella (2009) 
and Zhang and Chan (2009). The extension to the use of the full multivariate
affine GH distribution is detailed in Ghalanos et al (2011).
}
\keyword{classes}
