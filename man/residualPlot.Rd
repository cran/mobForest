\name{residualPlot}
\alias{residualPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Residual Diagnostics
}
\description{
Produces two plots: a) histogram of residuals, b) predicted Vs residuals.
This feature is applicable only when linear regression is considered as the node model. 
}
\usage{
residualPlot(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class 'mobForestOutput'
}
}
\details{
 Residuals are computed as difference between the predicted values of outcome (summarized across all trees) and observed values of outcome. The residual plots are typical when the fitted values are obtained through linear regression but not when logistic or Poisson regression is considered as a node model. Therefore, the residual plots are produced only when linear regression is considered. For logistic or Poisson models, a message is printed saying "Residual Plot not produced when logistic of Poisson regression is considered as the node model". 
}
\examples{
## Random Forest analysis of model based recursive partitioning

## load data  
data("BostonHousing", package = "mlbench")
\donttest{
## Recursive partitioning based on linear regression model medv ~ lstat with 3 trees. 
## 1 core/processor used. Supply more processors using 'processors' argument
rfout <- mobForestAnalysis(as.formula(medv ~ lstat), 
c("rad", "tax", "crim"), mobForest.controls = 
mobForest_control(ntree = 3, mtry = 2, replace = TRUE, 
alpha = 0.05, bonferroni = TRUE, minsplit = 120),
 data = BostonHousing, processors = 1, model = linearModel)

## Produces residual plot
residualPlot(rfout)
}
}
