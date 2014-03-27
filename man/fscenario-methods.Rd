\name{fscenario-methods}
\docType{methods}
\alias{fscenario}
\alias{fscenario-methods}
\alias{fscenario,ANY-method}
\title{Scenario Generation}
\description{
Generates a 1-ahead forecast scenario given a choice of data generating 
processes (for use in stochastic programming or risk management).
}
\usage{
fscenario(Data, sim = 1000, roll = 0, 
	model = c("gogarch", "dcc", "cgarch", "var", "mdist"), 
	spec = NULL, 
	var.model = list(lag = 1, lag.max = NULL, 
				lag.criterion = c("AIC", "HQ", "SC", "FPE"), 
				robust = FALSE, robust.control = list("gamma" = 0.25, 
				"delta" = 0.01, "nc" = 10, "ns" = 500)),
	mdist.model = list(distribution = c("mvn", "mvt", "manig"), 
				AR = TRUE, lag = 1),
	spd.control = list(lower = 0.1, upper = 0.9, type = "pwm", 
				kernel = "epanech"),
	cov.method = c("ML", "LW", "EWMA", "MVE", "MCD", "MVT", "BS"),
	cov.options = list(shrinkage=-1, lambda = 0.96),
	solver = "solnp", solver.control = list(), 
	fit.control = list(eval.se = FALSE), 
	cluster = NULL, save.output = FALSE, save.dir = getwd(),
	save.name = paste("S", sample(1:1000, 1), sep = ""), rseed  = NULL, ...)
}
\arguments{
\item{Data}{ An n-by-m data matrix or data.frame.}
\item{sim}{ The size of the simulated 1-ahead forecast.}
\item{roll}{ Whether to fit the data using (n - roll) periods and then
return a (roll+1) 1-ahead rolling simulated scenarios.}
\item{model}{ A choice of 5 models for generating scenarios.}
\item{spec}{ Required if choosing \sQuote{gogarch}, \sQuote{dcc} or \sQuote{cgarch}, 
in which case this represents a specification object (see rmgarch package) .}
\item{var.model}{ Required if model is var.}
\item{mdist.model}{ Required if model is mdist, and provides details for the model 
estimation (not yet implemented).}
\item{spd.control}{ Required if model is \dQuote{cgarch} and transformation is spd.}
\item{cov.method}{ For model \dQuote{var} this represents the choice of 
covariance matrix to use to generate random deviates.}
\item{cov.options}{ For model \dQuote{var} this provides the optional parameters
to certain types of covariance estimation methods.}
\item{solver}{ The choice of solver to use for all models but \dQuote{var},
and includes \sQuote{solnp}, \sQuote{nlminb} and \sQuote{nloptr}.}
\item{solver.control}{ Optional control options passed to the appropriate solver
chosen.}
\item{fit.control}{Control arguments passed to the fitting routine.}
\item{cluster}{ A cluster object created by calling \code{makeCluster} from 
the parallel package. If it is not NULL, then this will be used for parallel 
estimation of the refits (remember to stop the cluster on completion).}
\item{save.output}{ Whether output should be saved to file instead of being 
returned to the workspace.}
\item{save.dir}{ The directory to save output if save.output is TRUE.}
\item{save.name}{ The name of the file to save the output list.}
\item{rseed}{ A vector of length sim to initiate the random number generator.}
\item{...}{ Additional parameters passed to the model fitting routines. In 
particular, for the \sQuote{gogarch} model additional parameters are passed to
the ICA routines, whereas for the \sQuote{dcc} and \sQuote{cgarch} models 
this would include the \sQuote{realizedVol} xts matrix for the realGARCH model.}
}
\value{
A \code{\linkS4class{fScenario}} object containing the scenario and the model 
details (list). The scenario list contains a list of the
(roll+1) simulated forecast scenarios, the list of (roll+1) simulated forecast 
residuals, the forecast conditional mean, the forecast covariance and the list 
of random generator seed values used for replication. In addition, for the 
gogarch model the ICA whitening (K) and rotation matrices are also returned 
and required for replication of results (these may be entered in the 
\sQuote{gogarchspec} function). Use the \code{fitted} method on the object to
extract the simulated returns forecast.
}
\details{
The functionality here provides some wrapper functions, to create 1-ahead (and 
optionally rolling, useful for backtesting) scenarios for use in portfolio 
optimization using stochastic programming methods. The nature of these
data generating processes (as implemented here) and resulting optimization 
problems results in the so called anticipative class of stochastic programming 
models. If save.output is chosen, and given a save.dir, the scenario is saved 
(using save.name) and an object is returned containing an empty list for the 
scenario but with a model details list and the seed values. This can then be 
passed on to the \code{goload} function which can read from the directory and 
return a complete object with the scenario.
}
\author{Alexios Ghalanos}
\keyword{methods}