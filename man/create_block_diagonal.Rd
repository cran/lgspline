% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunctions.R
\name{create_block_diagonal}
\alias{create_block_diagonal}
\title{Create Block Diagonal Matrix}
\usage{
create_block_diagonal(matrix_list)
}
\arguments{
\item{matrix_list}{List of matrices to arrange diagonally}
}
\value{
Block diagonal matrix with input matrices on diagonal
}
\description{
Create Block Diagonal Matrix
}
\details{
Takes in a list of matrices, and returns a block-diagonal matrix with each
element of the list as one block. All off-diagonal elements are 0.
Matrices must have compatible dimensions.
}
\keyword{internal}
