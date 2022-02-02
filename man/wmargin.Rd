\name{wmargin}
\alias{wmargin}
\title{
Weighted Distribution Margin
}
\description{
Return the weighted margin of one of 3 elliptical distributions given a matrix 
of weights.
}
\usage{
wmargin(distribution = "mvnorm", weights, mean, Sigma, shape = NA, skew = NA)
}
\arguments{
  \item{distribution}{
  One of \sQuote{mvnorm}, \sQuote{mvlaplace} or \sQuote{mvt}.
}
  \item{weights}{
  Either a vector or matrix of weights, in the latter case must be of the same 
  row dimension as the covariance array.
}
  \item{mean}{
	Wither a vector or matrix of conditional distribution means, in the latter
	case must be of the same row dimension as the covariance array.
}
  \item{Sigma}{
  An array of covariances, usually returned by calling the \sQuote{rcov} method
  on one of the multivariate GARCH fitted objects.
}
  \item{shape}{
  The shape (d.o.f.) parameter of the multivariate student distribution.
}
 \item{skew}{
  Not currently required for the 3 distributions used.
 }
}
\details{
This is just a convenience function to return the weighted variance and mean of
the three elliptical distributions given a set of weights.
}
\value{
A matrix with each row representing the conditional weighted marginal density 
with corresponding parameters.
}
\author{Alexios Galanos}
