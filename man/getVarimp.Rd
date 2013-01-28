\name{getVarimp}
\alias{getVarimp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variable Importance
}
\description{
Variable importance scores computed through random forest analysis
}
\usage{
getVarimp(object)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{\linkS4class{mobForestOutput}} returned by \link[=mobForestAnalysis]{mobForestAnalysis()}
}
}
\references{
Leo Breiman (2001). Random Forests. \emph{Machine Learning}, 45(1), 5-32.\cr
}

\examples{
if(require("mlbench")) {

## recursive partitioning of a linear regression model

## load data  
data("BostonHousing", package = "mlbench")
BostonHousing <- BostonHousing[1:90,c("rad", "tax", "crim", "medv", "lstat")]

## Recursive partitioning based on linear regression model 
## medv ~ lstat with 3 trees. 
## 1 core/processor used. Supply more processors using 
## 'processors' argument
rfout <- mobForestAnalysis(as.formula(medv ~ lstat), 
c("rad", "tax", "crim"), mobForest.controls = 
mobForest_control(ntree = 3, mtry = 2, replace = TRUE, 
alpha = 0.05, bonferroni = TRUE, minsplit = 25), 
data = BostonHousing, processors = 1, model = linearModel)

## Returns a vector of variable importance scores 
getVarimp(rfout)		

}
}
