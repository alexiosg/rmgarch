\name{radical}
\alias{radical}
\title{
The Robust Accurate, Direct ICA aLgorithm (RADICAL).
}
\description{
An ICA algorithm based on an efficient entropy estimator (due to Vasicek) which 
is robust to outliers.
}
\usage{
radical(X, n.comp = dim(X)[2], demean = TRUE, pca.cov = c("ML", "LW", "ROB", "EWMA"), 
k = 150, augment = FALSE, replications = 30, sd = 0.175, firstEig = 1, 
lastEig = dim(X)[1], pcaE = NULL, pcaD = NULL, whiteSig = NULL, whiteMat = NULL, 
dewhiteMat = NULL, rseed = NULL, trace = FALSE, ...)
}
\arguments{
\item{X}{ The multidimensional signal matrix, where each column of matrix represents 
one observed signal.}
\item{n.comp}{ Number of independent components to estimate, defaults to the 
dimension of the data (rows). Is overwritten by \code{firstEig} and \code{lastEig}.}
\item{demean}{ (Logical) Whether the data should be centered.}
\item{pca.cov}{
The method to use for the calculation of the covariance matrix during the
PCA whitening phase. \dQuote{ML} is the standard maximum likelihood method,
\dQuote{LW} is the Ledoit-Wolf method, \dQuote{ROB} is the robust method from
the MASS package and \dQuote{EWMA} an exponentially weighted moving average
estimator. Optional parameters passed via the \dots argument.
}
\item{k}{ The number of angles at which to evaluate the contrast function. The 
ICA contrast function will be evaluated at K evenly spaced rotations from -Pi/4 
to Pi/4
}
\item{augment}{ Whether to augment the data (as explained in paper). For large 
datasets of >10,000 points this should be set to FALSE.}
\item{replications}{ This is the number of replicated points for each original 
point. The default value is 30. The larger the number of points in the data set, 
the smaller this value can be. For data sets of 10,000 points or more, point 
replication should be de-activated by setting augment to FALSE.}
\item{sd}{This is the standard deviation (noise) of the replicated points when 
using the augmentation option.}
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
\item{rseed}{Optionally provided seed to initialize the augmented data matrix.}
\item{trace}{To report progress in the console, set this to \code{TRUE}.}
\item{\dots}{Optional arguments passed to the pca.cov methods.}
}
\details{
The interested reader should consult the paper in the references section for 
details on the properties of the algorithm.\cr
The algorithm is quite slow, despite partial implementation in C++, and should 
only be used on small to medium sized sets.
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
Learned-Miller, A.G and Fisher III, J.W., 2003, ICA Using Spacings Estimates of 
Entropy, \emph{Journal of Machine Learning Research}, \bold{4}, 1271-1295.
\url{http://www.cs.umass.edu/~elm/ICA/}
}
\author{
Erik G. Learned-Miller for the Radical algorithm and Matlab package.\cr
Alexios Galanos for this R-port.
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
# unconditional covariance matrix) the data before applying the theICA algorithm.
IC <- radical(X, n.comp = 2)

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
