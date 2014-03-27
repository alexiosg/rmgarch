\name{goload-methods}
\docType{methods}
\alias{goload}
\alias{goload-methods}
\alias{goload,ANY-method}
\alias{goload,fScenario-method}
\alias{goload,fMoments-method}

\title{Load Scenario from File}
\description{
Loads a previously saved fScenario from file and returns
a \code{\linkS4class{fScenario}} or \code{\linkS4class{fMoments}} object.
}
\usage{
goload(object, ...)
}
\arguments{
\item{object}{ A \code{\linkS4class{fScenario}} or \code{\linkS4class{fMoments}} 
object which was created with save.output set to TRUE.}
\item{...}{ not used.}
}
\value{
A \code{\linkS4class{fScenario}} or \code{\linkS4class{fMoments}} object with 
the scenario or moments slot now filled with the saved data from file.
}
\details{
There are times when it is more efficient to save large scenarios to file 
(particularly when creating them in parallel), rather than returning them to the
user workspace. The save.output option in the \code{\link{fscenario}} and 
\code{\link{fmoments}} allows to do just that, returning instead a lighter 
object with an empty scenario slot, but with the model slot included, containing 
the details of the location and name of the saved scenario (or moments list). 
The goload function then takes this object, reads the location and name and 
loads the scenario (or moments) into its slot in the object and returns 
this to the user's workspace.
}
\author{Alexios Ghalanos}
\keyword{methods}