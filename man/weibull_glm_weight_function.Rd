% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunctions.R
\name{weibull_glm_weight_function}
\alias{weibull_glm_weight_function}
\title{Weibull GLM Weight Function for Constructing Information Matrix}
\usage{
weibull_glm_weight_function(
  mu,
  y,
  order_indices,
  family,
  dispersion,
  observation_weights,
  status
)
}
\arguments{
\item{mu}{Predicted survival times}

\item{y}{Observed response/survival times}

\item{order_indices}{Order of observations when partitioned to match "status" to "response"}

\item{family}{Weibull AFT family}

\item{dispersion}{Estimated dispersion parameter (\eqn{s^2})}

\item{observation_weights}{Weights of observations submitted to function}

\item{status}{Censoring indicator (1 = event, 0 = censored)
Indicates whether an event of interest occurred (1) or the observation was
right-censored (0). In survival analysis, right-censoring occurs when the
full survival time is unknown, typically because the study ended or the
subject was lost to follow-up before the event of interest occurred.}
}
\value{
Vector of weights for constructing the diagonal weight matrix \eqn{\textbf{W}}
in the information matrix \eqn{\textbf{G} = (\textbf{X}^{T}\textbf{W}\textbf{X} + \textbf{L})^{-1}}.
}
\description{
Computes diagonal weight matrix \eqn{\textbf{W}} for the information matrix
\eqn{\textbf{G} = (\textbf{X}^{T}\textbf{W}\textbf{X} + \textbf{L})^{-1}} in Weibull accelerated failure time (AFT) models.
}
\details{
This function generates weights used in constructing the information matrix
after unconstrained estimates have been found. Specifically, it is used in
the construction of the \eqn{\textbf{U}} and \eqn{\textbf{G}} matrices following initial unconstrained
parameter estimation.

These weights are analogous to the variance terms in generalized linear
models (GLMs). Like logistic regression uses \eqn{\mu(1-\mu)}, Poisson regression uses
\eqn{e^{\mu}}, and Linear regression uses constant weights, Weibull AFT models use
\eqn{\exp((\log y - \log \mu)/s)} where \eqn{s} is the scale (= \eqn{\sqrt{\text{dispersion}}}) parameter.
}
\examples{

## Demonstration of glm weight function in constrained model estimation
set.seed(1234)
n <- 1000
t1 <- rnorm(n)
t2 <- rbinom(n, 1, 0.5)

## Generate survival times
lambda <- exp(0.5 * t1 + 0.3 * t2)
yraw <- rexp(n, rate = 1/lambda)

## Introduce right-censoring
status <- rbinom(n, 1, 0.75)
y <- ifelse(status, yraw, runif(1, 0, yraw))

## Fit model demonstrating use of custom glm weight function
model_fit <- lgspline(y ~ spl(t1) + t2,
                      data.frame(y = y, t1 = t1, t2 = t2),
                      unconstrained_fit_fxn = unconstrained_fit_weibull,
                      family = weibull_family(),
                      need_dispersion_for_estimation = TRUE,
                      dispersion_function = weibull_dispersion_function,
                      glm_weight_function = weibull_glm_weight_function,
                      shur_correction_function = weibull_shur_correction,
                      status = status,
                      opt = FALSE,
                      K = 1)

print(summary(model_fit))


}
