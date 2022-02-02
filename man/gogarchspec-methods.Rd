\name{gogarchspec-methods}
\alias{gogarchspec}
\alias{gogarchspec-methods}
\alias{gogarchspec,ANY-method}
\title{function: GO-GARCH Specification}
\description{
Method for creating a GO-GARCH specification object prior to fitting.
}
\usage{
gogarchspec(mean.model = list(model = c("constant", "AR", "VAR"), robust = FALSE, 
lag = 1, lag.max = NULL, lag.criterion = c("AIC", "HQ", "SC", "FPE"), 
external.regressors = NULL, 
robust.control = list("gamma" = 0.25, "delta" = 0.01, "nc" = 10, "ns" = 500)), 
variance.model = list(model = "sGARCH", garchOrder = c(1,1), submodel = NULL, 
variance.targeting = FALSE), distribution.model = c("mvnorm", "manig", "magh"), 
ica = c("fastica", "radical"), 
ica.fix = list(A = NULL, K = NULL), ...) 
}
\arguments{
\item{mean.model}{
The mean specification. Allows for either a constant filtration of the return 
series, a univariate AR for each series with common lag (via the \dQuote{lag} 
argument) else a classical or robust Vector Autoregressive Model (VAR). The 
\sQuote{robust} option allows for a robust version of  VAR based on the 
multivariate Least Trimmed Squares Estimator  described in Croux and 
Joossens (2008). The \sQuote{robust.control} includes additional tuning 
parameters to the robust regression including the proportion to 
trim (\dQuote{gamma}), the critical value for Reweighted estimator (\dQuote{delta}), 
the number of subsets (\dQuote{ns}) and the number of C-steps (\dQuote{nc}).
The  external.regressors argument allows for a matrix of common external 
regressors in the constant, AR or VAR formulations.
}
\item{variance.model}{
The univariate variance specification for the independent factors of the 
GO-GARCH model.
}
  \item{distribution.model}{
The distributions supported are the multivariate normal (\dQuote{mvnorm}) and 
the multivariate affine NIG (\dQuote{manig}) and GHYP (\dQuote{magh}) 
distributions of Schmidt et al (see references).
}
\item{ica}{
The algorithm to use for extracting the independent components. 
The \code{\link{fastica}} and \code{\link{radical}} algorithms are the only
ICA algorithms currently allowed and locally implemented. See their 
documentation for a list of additional arguments possible, which may be passed
in the \code{\link{gogarchfit}} method.
}
\item{ica.fix}{
This allows the option of supplying the mixing matrix (A) and optionally the 
whitening Matrix (K). This is likely to be use when comparing different models 
(with the same mean filtration and dataset but different variance models) and 
you wish to use the same independent factors.
}
\item{\dots}{.}
}
\value{
A \code{\linkS4class{goGARCHspec}} object containing details of the GO-GARCH 
specification.
}
\author{Alexios Galanos}
\keyword{methods}
