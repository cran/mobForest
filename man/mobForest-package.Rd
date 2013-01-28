\name{mobForest-package}
\alias{mobForest-package}
\alias{mobForest}
\docType{package}
\title{
Random Forest methodology for model-based recursive partitioning

}
\description{
This package implements random forest method for model based recursive partitioning. 
}
\details{
\tabular{ll}{
Package: \tab mobForest\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2013-01-28\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
Depends: \tab parallel, party\cr
}

}
\author{
Nikhil Garge, Barry Eggleston, Georgiy Bobashev 
Maintainer: Nikhil Garge <ngarge@rti.org>
}
\references{
Achim Zeileis, Torsten Hothorn, and Kurt Hornik (2008). Model-Based Recursive Partitioning. \emph{Journal of Computational and Graphical Statistics}, 17(2), 492-514.\cr

Hothorn, T., Hornik, K. and Zeileis, A. (2006) Unbiased recursive partitioning: A conditional inference framework, \emph{J Comput Graph Stat}, 15, 651-674.\cr

Strobl, C., Malley, J. and Tutz, G. (2009) An introduction to recursive partitioning: rationale, application, and characteristics of classification and regression trees, bagging, and random forests, \emph{Psychol Methods}, 14, 323-348.\cr
}
\keyword{ package }
\examples{
\donttest{
set.seed(290875)
if(require("mlbench")) {
## Random Forest analysis of model based recursive partitioning

## load data  
data("BostonHousing", package = "mlbench")

## Recursive partitioning based on linear regression model medv ~ lstat with 3 trees. 
## 1 core/processor used. Supply more processors using 'processors' argument
rfout <- mobForestAnalysis(as.formula(medv ~ lstat), 
c("rad", "tax", "crim"), mobForest.controls = 
mobForest_control(ntree = 3, mtry = 2, replace = TRUE, 
alpha = 0.05, bonferroni = TRUE, minsplit = 120),
 data = BostonHousing, processors = 1, model = linearModel)

 rfout		## should produce descriptive output for mobForestOutput object
}
}
}
