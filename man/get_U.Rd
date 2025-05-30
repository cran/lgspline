% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunctions.R
\name{get_U}
\alias{get_U}
\title{Efficiently Construct U Matrix}
\usage{
get_U(G, A, K, nc, nca)
}
\arguments{
\item{G}{List of G matrices (\eqn{\textbf{G}})}

\item{A}{Constraint matrix (\eqn{\textbf{A}})}

\item{K}{Number of partitions minus 1 (\eqn{K})}

\item{nc}{Number of columns per partition}

\item{nca}{Number of constraint columns}
}
\value{
\eqn{\textbf{U}} matrix for constraints
}
\description{
Efficiently Construct U Matrix
}
\details{
Computes \eqn{\textbf{U} = \textbf{I} - \textbf{G}\textbf{A}(\textbf{A}^{T}\textbf{G}\textbf{A})^{-1}\textbf{A}^{T}} efficiently, avoiding unnecessary
multiplication of blocks of \eqn{\textbf{G}} with all-0 elements.
}
\keyword{internal}
