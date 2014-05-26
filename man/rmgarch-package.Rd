\name{rmgarch-package}
\alias{rmgarch-package}
\alias{rmgarch}
\title{The rmgarch package}
\description{
The rmgarch provides a selection of multivariate GARCH models with methods for
fitting, filtering, forecasting and simulation with additional support functions
for working with the returned objects. At present, the Generalized Orthogonal 
GARCH using Independent Components Analysis (\emph{ICA}) and Dynamic Conditional 
Correlation (with multivariate Normal, Laplace and Student distributions) models 
are fully implemented, with methods for spec, fit, filter, forecast, simulation, 
and rolling estimation and forecasting, as well as specialized functions to 
calculate and work with the weighted portfolio conditional density.
The Copula-GARCH model is also implemented with the multivariate Normal 
and Student distributions, with dynamic (DCC) and static estimation of 
the correlation.\cr
}
\details{
\tabular{ll}{
Package: \tab rmgarch\cr
Type: \tab Package\cr
Version: \tab 1.2-8\cr
Date: \tab 2014-05-27\cr
License: \tab GPL\cr
LazyLoad: \tab yes\cr
Depends: R (>= 3.0.2), methods, rugarch\cr
LinkingTo:\tab Rcpp, RcppArmadillo\cr
Imports: Rsolnp, MASS, parallel, Matrix, zoo, xts, Bessel, ff, fftw, shape, Kendall, spd, Rcpp\cr}
The main package functionality, currently supports the \emph{GO-GARCH} with ICA 
method, and is available through the \code{\link{gogarchspec}}, 
\code{\link{gogarchfit}}, \code{\link{gogarchfilter}}, \code{\link{gogarchforecast}}, 
\code{\link{gogarchsim}} and \code{\link{gogarchroll}} functions. The \emph{DCC} 
with multivariate Normal, Laplace and Student distributions is also supported 
with the main functionality in \code{\link{dccspec}}, \code{\link{dccfit}}, 
\code{\link{dccfilter}}, \code{\link{dccforecast}}, \code{\link{dccsim}} and
\code{\link{dccroll}}. The Normal and Student Copula-GARCH, with dynamic or 
static correlation, is implemented with the main functionality in 
\code{\link{cgarchspec}}, \code{\link{cgarchfit}}, \code{\link{cgarchfilter}}, 
and \code{\link{cgarchsim}}. Usual extractor and support methods for the 
multivariate GARCH models are documented in the class of the returned objects.\cr.
}
\section{How to cite this package}{
  Whenever using this package, please cite as\cr
 \preformatted{@Manual{Ghalanos_2014,
 author       = {Alexios Ghalanos},
 title        = {{rmgarch}: Multivariate GARCH models.},
 year         = {2014},
 note 	      = {R package version 1.2-8.}}}
}
\section{License}{
  The releases of this package is licensed under GPL version 3.
}
\author{Alexios Ghalanos}
\references{
Bollerslev, T. 1990, Modelling the coherence in short-run nominal exchange 
rates: a multivariate generalized ARCH model, \emph{The Review of Economics and 
Statistics}, \bold{72(3)}, 498--505.\cr
Broda, S.A. and Paolella, M.S. 2009, CHICAGO: A Fast and Accurate Method for 
Portfolio Risk Calculation, \emph{Journal of Financial Econometrics} \bold{7(4)}, 
412--436 .\cr
Cappiello, L., Engle, R.F. and Sheppard, K. 2006, Asymmetric dynamics in the 
correlations of global equity and bond returns, \emph{Journal of Financial 
Econometrics} \bold{4}, 537--572.\cr
Croux, C. and Joossens, K. 2008, Robust estimation of the vector autoregressive 
model by a least trimmed squares procedure, \emph{COMPSTAT}, 489--501.\cr
Chen, Y., Hardle, W., and Spokoiny, V. 2010, GHICA - Risk analysis with GH 
distributions and independent components, \emph{Journal of Empirical Finance},
\bold{17(2)}, 255--269.\cr
de Athayde, G.M. and Flores Jr, R.G. 2002, On Certain Geometric Aspects of 
Portfolio Optimisation with Higher Moments, \emph{mimeo}.\cr
Engle, R.F. 2002, Dynamic conditional correlation, \emph{Journal of Business and 
Economic Statistics} \bold{20}, 339--350.\cr
Engle, R.F. and Sheppard, K. 2001, Theoretical and empirical properties of 
dynamic conditional correlation multivariate GARCH, \emph{NBER Working Paper}.\cr
Genest, C., Ghoudi, K. and Rivest, L. 1995, A semiparametric estimation 
procedure of dependence parameters in multivariate families of distributions, 
\emph{Biometrika}, \bold{82}, 543--552.\cr
Ghalanos, A., Rossi, E., and Urga, G. (2013). \emph{Independent Factor 
Autoregressive Conditional Density Model}, \bold{forthcoming}.\cr
Paolella, M.S. 2007, Intermediate Probability - A Computational Approach, 
\emph{Wiley-Interscience}.\cr
Schmidt, R., Hrycej, T. and Stutzle 2006, Multivariate distribution models with 
generalized hyperbolic margins, \emph{Computational Statistics \& Data Analysis} 
\bold{50(8)}, 2065--2096.\cr
}