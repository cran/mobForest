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
data("BostonHousing", package = "mlbench")

rfout <- mobForestAnalysis(as.formula(medv ~ lstat + rm), c("zn", "indus", "chas", "nox", 
"age", "dis", "rad", "tax", "crim", "b", "ptratio"), mobForest.controls = 
mobForest_control(ntree = 8, mtry = 4, replace = TRUE, fraction = 0.632, minsplit = 80),
data = BostonHousing, 
processors = 1,
model = linearModel)

## prints predictive accuracy output 
pacc <- PredictiveAccuracy(rfout)
pacc

}
