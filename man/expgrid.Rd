% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunctions.R
\name{expgrid}
\alias{expgrid}
\title{Generate Grid Indices Without expand.grid()}
\usage{
expgrid(vec_list, indices)
}
\arguments{
\item{vec_list}{List of vectors to combine}

\item{indices}{Indices of combinations to return}
}
\value{
Data frame of selected combinations
}
\description{
Generate Grid Indices Without expand.grid()
}
\details{
Returns selected combinations from the cartesian product of \code{vec_list}
without constructing full \code{expand.grid()} for memory efficiency.
}
\keyword{internal}
