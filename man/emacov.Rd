\name{ewmacov}
\alias{ewmacov}
\title{
Exponentially Weighted Moving Average (EWMA) Covariance Filter
}
\description{
A fast EWMA covariance and mean filter, including diagonal variants.
}
\usage{
ewmacov(X, lambda.v = 0.96, lambda.m = 0.96, cor = FALSE)
}
\arguments{
\item{X}{ A multidimensional xts matrix.}
\item{lambda.v}{Either a scalar or vector of covariance decay coefficients.}
\item{lambda.m}{ Either a scalar or vector of mean decay coefficients. If NULL,
then it will be assumed that X is already centered and no additional centering
will take place.}
\item{cor}{(logical) Whether to return the correlation matrix (in addition to
the covariance matrix)}.
}
\details{
}
\value{
A list containing the following values:
\item{cov}{The filtered array of time varying covariances.}
\item{cor}{The filtered array of time varying correlations (if cor was TRUE).}
\item{M}{The filtered xts matrix of the (time varying) means 
(if lambda.m was not NULL).}
}
\author{
Alexios Ghalanos
}
\examples{
\dontrun{
}}
\keyword{multivariate}