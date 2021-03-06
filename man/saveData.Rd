% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methylInheritanceSimInternalMethods.R
\name{saveData}
\alias{saveData}
\title{Save data created during the simulation}
\usage{
saveData(pathOut, extension, gRanges, methylData, methUnit, diffData,
  saveGRanges, saveMethylKit, runAnalysis)
}
\arguments{
\item{pathOut}{a string of \code{character}, the path 
where the files are saved.}

\item{extension}{a string of \code{character} representing the extension
that will be given to the saved files.}

\item{gRanges}{a \code{list} of \code{methylRawList}  TODO}

\item{methylData}{a \code{list} of \code{methylRawList}, the results
of the normalization of the coverage.}

\item{methUnit}{a \code{list} of \code{methylBase}, the results of the
base filtering for all samples.}

\item{diffData}{a \code{list} of \code{methylDiff}, the results of the
calculation of differential methylation statistics.}

\item{saveGRanges}{a \code{logical}, when \code{true}, files containing 
\code{GRangeaList} are saved.}

\item{saveMethylKit}{a \code{logical}, when \code{TRUE}, files 
\code{methylRawList} object are saved.}

\item{runAnalysis}{a \code{logical}, when \code{TRUE}, two files related
to the analysis are saved.}
}
\value{
\code{0} indicating that the function has been successful.
}
\description{
Save data created during the simulation.
}
\examples{

## TODO

}
\author{
Pascal Belleau, Astrid Deschenes
}
\keyword{internal}
