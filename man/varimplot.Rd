\name{varimplot}
\alias{varimplot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variable importance plot
}
\description{
A plot with variable importance score on X-axis and variable name on Y-axis.
}
\usage{
varimplot(object)
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

\seealso{
\link[=getVarimp]{getVarimp}
}
\examples{
\donttest{
if(require("mlbench")) {

## recursive partitioning of a linear regression model
## load data  
data("BostonHousing", package = "mlbench")

## Recursive partitioning based on linear regression model medv ~ lstat with 3 trees. 
## 1 core/processor used. Supply more processors using 'processors' argument
rfout <- mobForestAnalysis(as.formula(medv ~ lstat), 
c("rad", "tax", "crim"), mobForest.controls = 
mobForest_control(ntree = 3, mtry = 2, replace = TRUE, 
alpha = 0.05, bonferroni = TRUE, minsplit = 120),
 data = BostonHousing, processors = 1, model = linearModel)

## Produce variable importance plot
varimplot(rfout)  
}
}
}
