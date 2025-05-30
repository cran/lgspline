% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunctions.R
\name{unconstrained_fit_weibull}
\alias{unconstrained_fit_weibull}
\title{Unconstrained Weibull Accelerated Failure Time Model Estimation}
\usage{
unconstrained_fit_weibull(
  X,
  y,
  LambdaHalf,
  Lambda,
  keep_weighted_Lambda,
  family,
  tol = 1e-08,
  K,
  parallel,
  cl,
  chunk_size,
  num_chunks,
  rem_chunks,
  order_indices,
  weights,
  status
)
}
\arguments{
\item{X}{Design matrix of predictors}

\item{y}{Survival/response times}

\item{LambdaHalf}{Square root of penalty matrix (\eqn{\boldsymbol{\Lambda}^{1/2}})}

\item{Lambda}{Penalty matrix (\eqn{\boldsymbol{\Lambda}})}

\item{keep_weighted_Lambda}{Flag to retain weighted penalties}

\item{family}{Distribution family specification}

\item{tol}{Convergence tolerance (default 1e-8)}

\item{K}{Number of partitions minus one (\eqn{K})}

\item{parallel}{Flag for parallel processing}

\item{cl}{Cluster object for parallel computation}

\item{chunk_size}{Processing chunk size}

\item{num_chunks}{Number of computational chunks}

\item{rem_chunks}{Remaining chunks}

\item{order_indices}{Observation ordering indices}

\item{weights}{Optional observation weights}

\item{status}{Censoring status indicator (1 = event, 0 = censored)
Indicates whether an event of interest occurred (1) or the observation was
right-censored (0). In survival analysis, right-censoring occurs when the
full survival time is unknown, typically because the study ended or the
subject was lost to follow-up before the event of interest occurred.}
}
\value{
Optimized beta parameter estimates (\eqn{\boldsymbol{\beta}}) for Weibull AFT model
}
\description{
Estimates parameters for an unconstrained Weibull accelerated failure time
(AFT) model supporting right-censored survival data.

This both provides a tool for actually fitting Weibull AFT Models, and
boilerplate code for users who wish to incorporate Lagrangian multiplier
smoothing splines into their own custom models.
}
\details{
Estimation Approach:
The function employs a two-stage optimization strategy for fitting
accelerated failure time models via maximum likelihood:

1. Outer Loop: Estimate Scale Parameter using Brent's method

2. Inner Loop: Estimate Regression Coefficients (given scale) using
   damped Newton-Raphson.
}
\examples{

## Simulate survival data with covariates
set.seed(1234)
n <- 1000
t1 <- rnorm(n)
t2 <- rbinom(n, 1, 0.5)

## Generate survival times with Weibull-like structure
lambda <- exp(0.5 * t1 + 0.3 * t2)
yraw <- rexp(n, rate = 1/lambda)

## Introduce right-censoring
status <- rbinom(n, 1, 0.75)
y <- ifelse(status, yraw, runif(1, 0, yraw))
df <- data.frame(y = y, t1 = t1, t2 = t2)

## Fit model using lgspline with Weibull AFT unconstrained estimation
model_fit <- lgspline(y ~ spl(t1) + t2,
                      df,
                      unconstrained_fit_fxn = unconstrained_fit_weibull,
                      family = weibull_family(),
                      need_dispersion_for_estimation = TRUE,
                      dispersion_function = weibull_dispersion_function,
                      glm_weight_function = weibull_glm_weight_function,
                      shur_correction_function = weibull_shur_correction,
                      status = status,
                      opt = FALSE,
                      K = 1)

## Print model summary
summary(model_fit)

}
\keyword{internal}
