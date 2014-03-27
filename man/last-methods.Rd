\name{last-methods}
\docType{methods}
\alias{last-methods}
\alias{last,ANY-method}
\alias{last,array-method}
\alias{first-methods}
\alias{first,ANY-method}
\alias{first,array-method}
\alias{first}
\alias{last}
\title{First and Last methods for accessing objects}
\description{
Functions for accessing first-n and last-n values of an object 
(similar to head and tail).
}
\usage{
last(x, index = 1, ...)
first(x, index = 1, ...)
}
\arguments{
\item{x}{ Currently only arrays supported.}
\item{index}{ First or Last n-indices to return values for.}
\item{...}{ For expansion to other classes.}
}
\section{Methods}{
\describe{
\item{\code{signature(x = "ANY")}}{}
\item{\code{signature(x = "array")}}{}
}}
\author{Alexios Ghalanos}
\keyword{methods}