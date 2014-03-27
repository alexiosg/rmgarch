\name{goGARCHfft-class}
\docType{class}
\alias{goGARCHfft-class}
\alias{dfft,goGARCHfft-method}
\alias{dfft}
\alias{pfft,goGARCHfft-method}
\alias{pfft}
\alias{qfft,goGARCHfft-method}
\alias{qfft}
\alias{nportmoments,goGARCHfft-method}
\alias{nportmoments}
\title{Class: GO-GARCH portfolio density}
\description{Class for the GO-GARCH portfolio density}
\section{Objects from the Class}{
The class is returned by calling the function \code{\link{convolution}}  on 
objects of class \code{\linkS4class{goGARCHfit}}, \code{\linkS4class{goGARCHfilter}}, 
\code{\linkS4class{goGARCHforecast}}, \code{\linkS4class{goGARCHsim}} and
\code{\linkS4class{goGARCHroll}}
}
\section{Slots}{
  \describe{
    \item{\code{dist}:}{ A list with the portfolio density and other details.}
    \item{\code{model}:}{ A list with the model details carried across objects.}
  }
}
\section{Methods}{
\describe{
\item{dfft}{\code{signature(object = "goGARCHfft")}:
The takes additional argument \dQuote{index} to indicate the particular time 
point, and returns an interpolated density function which may be called like 
any other \dQuote{d} type density function.
}
\item{pfft}{\code{signature(object = "goGARCHfft")}
The takes additional argument \dQuote{index} to indicate the particular time 
point, and returns an interpolated distribution function which may be called 
like any other \dQuote{p} type distribution function.}
\item{qfft}{\code{signature(object = "goGARCHfft")} 
This takes additional argument \dQuote{index} to indicate the particular time 
point, and returns an interpolated quantile function which may be called like 
any other \dQuote{q} type quantile function. This may also be used to generate 
pseudo-random variables from the distribution by using random standard uniform 
numbers as inputs.}
\item{nportmoments}{\code{signature(object = "goGARCHfft")}:
Calculate and returns a matrix of the first 4 standardized moments by evaluation 
of the portfolio density using quadrature based method (i.e. calling R's 
\dQuote{integrate} function on the portfolio FFT based density). Depending on
the GOGARCH class the density was based (e.g. goGARCHfit vs goGARCHforecast), 
the format of the output will be different, and generally follow the 
format \sQuote{rules} of that class.}
}
}
\section{notes}{
In the case that \code{convolution} was called on a \code{\linkS4class{goGARCHforecast}}
or \code{\linkS4class{goGARCHroll}} object, the \code{dist} slot will contain the max 
of n.ahead or n.roll. There should be no confusion here since the multivariate 
forecast methods in rmgarch only allow either n.ahead>1 with n.roll = 0 
(pure unconditional), or n.ahead = 1 with n.roll>=0 (pure rolling), and only the
latter in the case of a \code{\link{gogarchroll}}. While the \code{nportmoments} 
method reconstitutes the forecasts into a more familiar form (n.ahead x n.moments x 
(n.roll+1)), this does not make sense for the distribution methods (d*, p*, 
and q*), and it is understood that when the user calls for example 
\code{dfft(object, index=5)} on an object created from a forecast with 
n.ahead=10 and n.roll=0, the index is meant to indicate the unconditional 
density forecast at time T+5. Similarly, when calling 
code{dfft(object, index=0)} on an object created from a forecast with 
n.ahead=1 and n.roll = 1 (remember that n.roll is zero based), the index is 
meant to indicate the first (of two, since rolls = 0:1) rolling forecast 
density.
}
\author{Alexios Ghalanos}
\keyword{classes}