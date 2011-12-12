\name{prediction_output}
\alias{prediction_output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions and predictive accuracy
}
\description{
Predictions and predictive accuracy estimates
}
\usage{
prediction_output(predMean = numeric(), predSd = numeric(), 
residual = numeric(), R2 = numeric(), mse = numeric(), 
overallR2 = numeric(), predType = character())
}
\arguments{
  \item{predMean}{
Mean predictions across trees.
}
  \item{predSd}{
Standard deviation predictions across trees.
}
  \item{residual}{
Residuals (predicted outcome - observed outcome).
}
  \item{R2}{
Predictive accuracy across trees
}
  \item{mse}{
MSE across trees
}
  \item{overallR2}{
Overall R2
}
  \item{predType}{
Out-of-bag data or test data or learning data.
}
}
\details{
This function takes predictions and predictive accuracy estimates as input arguments and creates objects of class \code{\linkS4class{predictionOutput}}.
}
\value{
An object of class \code{\linkS4class{predictionOutput}}.
}

\seealso{
\code{\linkS4class{predictionOutput}}, \code{\link{mobForestAnalysis}}
}
