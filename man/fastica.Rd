\name{fastica}
\alias{fastica}
\title{
Fast Fixed Point ICA
}
\description{
The fast fixed point algorithm for independent component analysis and projection 
pursuit based on the direct translation to R of the FastICA program of the 
original authors at the Helsinki University of Technology.
}
\usage{
fastica(X, approach = c("symmetric", "deflation"), n.comp = dim(X)[2], demean = TRUE, 
pca.cov = c("ML", "LW", "ROB", "EWMA"), gfun = c("pow3", "tanh", "gauss", "skew"), 
finetune = c("none", "pow3", "tanh", "gauss", "skew"), tanh.par = 1, gauss.par = 1, 
step.size = 1, stabilization = FALSE, epsilon = 1e-4, maxiter1 = 1000, maxiter2 = 5, 
A.init = NULL, pct.sample = 1, firstEig = NULL, lastEig = NULL, 
pcaE = NULL, pcaD = NULL, whiteSig = NULL, whiteMat = NULL, dewhiteMat = NULL, 
rseed = NULL, trace = FALSE, ...)
}
\arguments{
\item{X}{ The multidimensional signal matrix, where each column of matrix 
represents one observed signal.}
\item{approach}{ The decorrelation approach to use, with \dQuote{symmetric} 
estimating the components in parallel while \dQuote{deflation} estimating 
one-by-one as in projection pursuit.
}
\item{n.comp}{ Number of independent components to estimate, defaults to the 
dimension of the data (rows). Is overwritten by \code{firstEig} and \code{lastEig}.}
\item{demean}{ (Logical) Whether the data should be centered.}
\item{pca.cov}{ The method to use for the calculation of the covariance matrix 
during the PCA whitening phase. \dQuote{ML} is the standard maximum likelihood 
method, \dQuote{LW} is the Ledoit and Wolf method, \dQuote{ROB} is the robust method 
from the MASS package and \dQuote{EWMA} an exponentially weighted moving average
estimator. Optional parameters passed via the (...) argument.
}
\item{gfun}{ The nonlinearity algorithm to use in the fixed-point algorithm.}
\item{finetune}{ The nonlinearity algorithm for fine-tuning.}
\item{tanh.par}{ Control parameter used when nonlinearity algorithm 
equals \dQuote{tanh}.}
\item{gauss.par}{ Control parameter used when nonlinearity algorithm equals 
\dQuote{gauss}.}
\item{step.size}{ Step size. If this is anything other than 1, the program will 
use the stabilized version of the algorithm.}
\item{stabilization}{ Controls whether the program uses the stabilized version 
of the algorithm. If the stabilization is on, then the value of \code{step.size} 
can momentarily be halved if the program estimates that the algorithm is stuck 
between two points (this is called a stroke). Also if there is no convergence 
before half of the maximum number of iterations has been reached then the 
\code{step.size} will be halved for the rest of the rounds.
}
\item{epsilon}{ Stopping criterion. Default is 0.0001.}
\item{maxiter1}{ Maximum number of iterations for \code{gfun} algorithm.}
\item{maxiter2}{ Maximum number of iterations for \code{finetune} algorithm.}
\item{A.init}{ Initial guess for the mixing matrix A. Defaults to a random 
(standard normal) filled matrix (no.signals by no.factors).}
\item{pct.sample}{ Percentage [0-1] of samples used in one iteration. Samples are 
chosen at random.}
\item{firstEig}{ This and \code{lastEig} specify the range for eigenvalues 
that are retained, \code{firstEig} is the index of largest eigenvalue to be 
retained. Making use of this option overwrites \code{n.comp}.
}
\item{lastEig}{ This is the index of the last (smallest) eigenvalue to be 
retained and overwrites \code{n.comp} argument.}
\item{pcaE}{ Optionally provided eigenvector (must also supply \code{pcaD}).}
\item{pcaD}{ Optionally provided eigenvalues (must also supply \code{pcaE}).}
\item{whiteSig}{Optionally provided Whitened signal.}
\item{whiteMat}{ Optionally provided Whitening matrix (no.factors by no.signals).}
\item{dewhiteMat}{Optionally provided dewhitening matrix (no.signals by no.factors).}
\item{rseed}{ Optionally provided seed to initialize the mixing matrix A 
(when \code{A.init} not provided).}
\item{trace}{To report progress in the console, set this to \sQuote{TRUE}.}
\item{\dots}{Optional arguments passed to the pca.cov methods.}
}
\details{
The fastica program is a direct translation into R of the FastICA Matlab program 
of Gaevert, Hurri, Saerelae, and Hyvaerinen with some extra features. All 
computations are currently implemented in R so for very large dimensional sets
alternative implementations may be faster. Porting part of the code to C++ may be
implemented in a future version.
}
\value{
A list containing the following values:
\item{A}{Estimated Mixing Matrix (no.signals by no.factors).}
\item{W}{Estimated UnMixing Matrix (no.factors by no.signals). }
\item{U}{Estimated rotation Matrix (no.factors by no.factors). }
\item{S}{The column vectors of estimated independent components (no.obs by no.factors).}
\item{C}{Estimated Covariance Matrix (no.signals by no.signals).}
\item{whiteningMatrix}{The Whitening matrix (no.factors by no.signals).}
\item{dewhiteningMatrix}{The de-Whitening matrix  (no.signals by no.factors).}
\item{rseed}{The random seed used (if any) for initializing the mixing matrix A.}
\item{elapsed}{The elapsed time.}
}
\references{
Hyvaerinen, A. and Oja,.E , 1997, A fast fixed-point algorithm for independent 
component analysis, \emph{Neural Computation}, \bold{9(7)}, 1483-1492. Reprinted in 
\emph{Unsupervised Learning}, G. Hinton and T. J. Sejnowski, 1999, MIT Press.\cr
}
\note{
Since version 1.0-3 the multidimensional signal matrix is now the usual row by 
column matrix, where the rows represent observations and columns the signals. 
Before this version, the reverse was true in keeping with the original version 
of the program.\cr
Dimensionality reduction can be achieved in the PCA stage by use of either 
\code{n.comp} in which case the n.comp largest eigenvalues are chosen, else by
selection of \code{firstEig} and \code{lastEig} which overwrites the choice of
\code{n.comp}.
}
\author{
Hugo Gaevert, Jarmo Hurri, Jaakko Saerelae, and Aapo Hyvaerinen for the original 
FastICA package for matlab.\cr 
Alexios Ghalanos for this R-port.
}
\examples{
\dontrun{
# create a set of independent signals S, glued together by a mixing matrix A
# (note the notation and matrix multiplication direction as we are dealing with
# row rather than column vectors)
set.seed(100)
S <- matrix(runif(10000), 5000, 2)
A <- matrix(c(1, 1, -1, 2), 2, 2, byrow = TRUE)
# the mixed signal X
X = S \%*\% t(A)
# The function centers and whitens (by the eigenvalue decomposition of the 
# unconditional covariance matrix)  the data before applying the theICA algorithm.
IC <- fastica(X, n.comp = 2, approach = "symmetric", gfun = "tanh", trace  = TRUE, 
A.init = diag(2))

# demeaned data:
X_bar = scale(X, scale = FALSE)

# whitened data:
X_white = X_bar \%*\% t(IC$whiteningMatrix)

# check whitening:
# check correlations are zero
cor(X_white)
# check diagonals are 1 in covariance
cov(X_white)

# check that the estimated signals(S) multiplied by the
# estimated mxing matrix (A) are the same as the original dataset (X)
round(head(IC$S \%*\% t(IC$A)), 12) == round(head(X), 12)

# do some plots:
par(mfrow = c(1, 3))
plot(IC$S \%*\% t(IC$A), main = "Pre-processed data")
plot(X_white, main = "Whitened and Centered components")
plot(IC$S, main = "ICA components")
}}
\keyword{multivariate}