% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HelperFunctions.R
\name{weibull_family}
\alias{weibull_family}
\title{Weibull Family for Survival Model Specification}
\usage{
weibull_family()
}
\value{
A list containing family-specific components for survival model estimation
}
\description{
Creates a compatible family object for Weibull accelerated failure time (AFT)
models with customizable tuning options.

This both provides a tool for actually fitting Weibull AFT Models, and
boilerplate code for users who wish to incorporate Lagrangian multiplier
smoothing splines into their own custom models.
}
\details{
Provides a comprehensive family specification for Weibull AFT models, including Family
name, link function, inverse link function, and custom loss function for model tuning

Supports right-censored survival data with flexible parameter estimation.
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

## Prepare data
df <- data.frame(y = y, t1 = t1, t2 = t2, status = status)

## Fit model using custom Weibull family
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

summary(model_fit)

}
