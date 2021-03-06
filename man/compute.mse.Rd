% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute.mse.R
\name{compute.mse}
\alias{compute.mse}
\title{Predictive accuracy estimates (MSE) across trees for linear or poisson
regression model.}
\usage{
compute.mse(response, predictions)
}
\arguments{
\item{response}{A vector of actual response of outcome variable.}

\item{predictions}{A vector of predicted response for the same outcome
variable.}
}
\value{
MSE estimates
}
\description{
Predictive accuracy estimates (MSE) across trees for linear or poisson
regression model.
}
\examples{
# The MSE should be 2.5. Off by 2 half the time, off by 1 the other half
response <- matrix(c(rep(0,100), rep(10,100)))
predictions <-
    matrix(nrow=20, ncol = 3,
           data = c(rep(1,100), rep(8,100), rep(1,100), rep(8,100),
                    rep(1,100), rep(8,100)))
compute.mse(response, predictions)

}
