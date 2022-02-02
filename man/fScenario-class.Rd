\name{fScenario-class}
\docType{class}
\alias{fScenario-class}
\alias{show,fScenario-method}
\alias{goget,fScenario-method}
\alias{goget,ANY-method}
\alias{goget}
\alias{fitted,fScenario-method}
\title{Class \code{"fScenario"}}
\description{
Object returned from calling \code{\link{fscenario}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("fScenario", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{scenario}:}{Object of class \code{"vector"} A list with the
    (roll+1) scenario matrices.}
    \item{\code{model}:}{Object of class \code{"vector"} A list with details of 
    data generating process.}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "fScenario")}: Summary method. }
    \item{goget}{\code{signature(object = "fScenario")}:  Get a specified 
    \sQuote{arg} from the object (only \sQuote{scenario} used).}
    \item{fitted}{\code{signature(object = "fScenario")}: Returns an array
    of the simulated scenario returns, of dimensions n.sim by n.assets 
    by (roll+1), with third dimension labels the actual forecast time index,
    and second dimension labels the asset names. The last forecast scenario
    will always be completely out of sample so the time index label for that 
    is generated using the \code{generatefwd} function in the rugarch package.}
    }
}
\author{Alexios Galanos}
\keyword{classes}
