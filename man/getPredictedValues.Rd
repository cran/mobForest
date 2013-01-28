\name{getPredictedValues}
\alias{getPredictedValues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions
}
\description{
Get predictions summarized across trees for out-of-bag cases or all cases or cases from new test data.
}
\usage{
getPredictedValues(object, OOB = TRUE, newdata = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{mobForestOutput}
}
  \item{OOB}{
a logical determining whether to return predictions from the out-of-bag sample or the learning sample (not suggested).
}
  \item{newdata}{
a logical determining whether to return predictions from test data. If newdata = TRUE, then OOB argument is ignored.
}
}
\value{
matrix with three columns: 1) Mean Predictions across trees, 2)  Standard deviation of predcitions across trees, and 3) Residual (mean predicted - observed). The third column is applicable only when linear regression is considered as the node model.
}
\examples{
## Random Forest analysis of model based recursive partitioning

## load data  
data("BostonHousing", package = "mlbench")
BostonHousing <- BostonHousing[1:90,c("rad", "tax", "crim", "medv", "lstat")]

## Recursive partitioning based on linear regression model medv ~ lstat with 3 trees. 
## 1 core/processor used. Supply more processors using 'processors' argument
rfout <- mobForestAnalysis(as.formula(medv ~ lstat), 
c("rad", "tax", "crim"), mobForest.controls = 
mobForest_control(ntree = 3, mtry = 2, replace = TRUE, 
alpha = 0.05, bonferroni = TRUE, minsplit = 25), 
data = BostonHousing, processors = 1, model = linearModel)

## Obtain out-of-bag predicted values
OOB.predMat <- getPredictedValues(rfout, OOB = TRUE)
OOB.pred = OOB.predMat[,1]
}
