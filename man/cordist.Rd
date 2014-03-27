\name{cordist}
\alias{cordist}
\title{
A Correlation Distance Measure
}
\description{
Given an array of correlation snapshots in time, returns a matrix of some rolling
distance measure on the correlations.
}
\usage{
cordist(R, distance = c("ma","ms","meda","meds","eigen", "cmd"), n = 25, 
plot = TRUE, dates = NULL, title = NULL)
}
\arguments{
  \item{R}{
  An array of correlations.
}
  \item{distance}{
  The measure to use to capture the change between 2 correlation matrices 
  (see details).
}
  \item{n}{
  The distance between 2 correlation matrices.
}
  \item{plot}{
  Whether to create a heatmap plot of the result.
}
  \item{dates}{
   A \code{POSIXct} vector of dates to use for the heatmap (recommend to
   supply).
}
 \item{title}{
 Title for the heatmap plot.
 }
}
\details{
This function provides for a visualization of dynamic correlation distance 
between periods with a number of plausible measures including \dQuote{ma} 
(mean absolute), \dQuote{ms} (mean squared), \dQuote{meda} (median absolute), 
\dQuote{meds} (median squared) \dQuote{eigen} (largest eigenvalue difference) and
\dQuote{cmd} (correlation matrix distance). See the references for more details. 
}
\value{
A symmetric matrix of the rolling distance measure for each period.
}
\references{
Munnix, M. C., Shimada, T., Schafer, R., Leyvraz, F., Seligman, T. H., 
Guhr, T., & Stanley, H. E. (2012). Identifying states of a financial market. 
\emph{Scientific Reports 2}.\cr
Herdin, M., Czink, N., Ozcelik, H., & Bonek, E. (2005). Correlation matrix 
distance, a meaningful measure for evaluation of non-stationary MIMO channels. 
\emph{Vehicular Technology Conference, 2005, IEEE 61st}, \bold{1}, 136--140.
}
\author{Alexios Ghalanos}