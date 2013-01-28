\name{print.predAccuracyEstimates}
\alias{print.predAccuracyEstimates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictive Accuracy Report
}
\description{
Summarizes predictive Accuracy
}
\usage{
\method{print}{predAccuracyEstimates}(x, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class 'predAccuracyEstimates' return by \code{"\link[=PredictiveAccuracy]{PredictiveAccuracy()}"} function
}
  \item{\dots}{
Additional arguments to print method
}
}
\examples{
\donttest{
data("BostonHousing", package = "mlbench")

rfout <- mobForestAnalysis(as.formula(medv ~ lstat), 
c("rad", "tax", "crim"), mobForest.controls = 
mobForest_control(ntree = 3, mtry = 2, replace = TRUE, 
alpha = 0.05, bonferroni = T, minsplit = 120), 
data = BostonHousing, processors = 1, model = linearModel)

## prints predictive accuracy output 
pacc <- PredictiveAccuracy(rfout)
pacc
}
}
