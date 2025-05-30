% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunctions.R
\name{matsqrt}
\alias{matsqrt}
\title{Calculate Matrix Square Root}
\usage{
matsqrt(mat)
}
\arguments{
\item{mat}{A symmetric, positive-definite matrix \eqn{\textbf{M}}}
}
\value{
A matrix \eqn{\textbf{B}} such that \eqn{\textbf{B}\textbf{B} = \textbf{M}}
}
\description{
Calculate Matrix Square Root
}
\details{
For matrix \eqn{\textbf{M}}, computes \eqn{\textbf{B}} where \eqn{\textbf{B}\textbf{B} = \textbf{M}} using eigenvalue decomposition:

1. Compute eigendecomposition \eqn{\textbf{M} = \textbf{V}\textbf{D}\textbf{V}^T}

2. Set eigenvalues below \code{sqrt(.Machine$double.eps)} to 0 for stability

3. Take elementwise square root of eigenvalues: \eqn{\textbf{D}^{1/2}}

4. Reconstruct as \eqn{\textbf{B} = \textbf{V} \textbf{D}^{1/2} \textbf{V}^T}

This provides the unique symmetric positive-definite square root.

You can use this to help construct a custom \code{Vhalf_fxn} for fitting
correlation structures, see \code{\link{lgspline}}.
}
\examples{
## Identity matrix
m1 <- diag(2)
matsqrt(m1)  # Returns identity matrix

## Compound symmetry correlation matrix
rho <- 0.5
m2 <- matrix(rho, 3, 3) + diag(1-rho, 3)
B <- matsqrt(m2)
# Verify: B \%**\% B approximately equals m2
all.equal(B \%**\% B, m2)

## Example for correlation structure
n_blocks <- 2  # Number of subjects
block_size <- 3  # Measurements per subject
rho <- 0.7  # Within-subject correlation
# Correlation matrix for one subject
R <- matrix(rho, block_size, block_size) +
     diag(1-rho, block_size)
# Full correlation matrix for all subjects
V <- kronecker(diag(n_blocks), R)
Vhalf <- matsqrt(V)

}
