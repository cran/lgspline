% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunctions.R
\name{AGAmult_wrapper}
\alias{AGAmult_wrapper}
\title{Efficient Matrix Multiplication for \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}}
\usage{
AGAmult_wrapper(
  G,
  A,
  K,
  nc,
  nca,
  parallel,
  cl,
  chunk_size,
  num_chunks,
  rem_chunks
)
}
\arguments{
\item{G}{List of G matrices (\eqn{\textbf{G}})}

\item{A}{Constraint matrix (\eqn{\textbf{A}})}

\item{K}{Number of partitions minus 1 (\eqn{K})}

\item{nc}{Number of columns per partition}

\item{nca}{Number of constraint columns}

\item{parallel}{Use parallel processing}

\item{cl}{Cluster object}

\item{chunk_size}{Chunk size for parallel}

\item{num_chunks}{Number of chunks}

\item{rem_chunks}{Remaining chunks}
}
\value{
Matrix product \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}
}
\description{
Efficient Matrix Multiplication for \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}}
}
\details{
Computes \eqn{\textbf{A}^{T}\textbf{G}\textbf{A}} efficiently in parallel chunks using \code{AGAmult_chunk()}.
}
\keyword{internal}
