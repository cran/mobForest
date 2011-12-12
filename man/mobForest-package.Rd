\name{mobForest-package}
\alias{mobForest-package}
\alias{mobForest}
\docType{package}
\title{
Random Forest implementation of model-based recursive partitioning
}
\description{
'mobForest' package implements random forest method of recursive partitioning based on  parametric models. 'mobForest' uses the code in party package, developed by Zeileis et al (2008). The function mob(), within party package, is modified to implement the random forest version of model based recursive partitioning. The main input function mobForestAnalysis() takes all input parameters to construct trees, compute out-of-bag errors, predictions, and overall accuracy of forest. The algorithm performs parallel computation using clusterApply() function within snow package. The trees are constructed across the nodes/cores simultaneously, improving the efficiency of the algorithm.
}
\details{
\tabular{ll}{
Package: \tab mobForest\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-12-12\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
Depends: \tab parallel, party\cr
}

}
\author{
Nikhil Garge, Barry Eggleston 
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
c("rad", "tax", "crim"),
mobForest.controls = 
mobForest_control(ntree = 3, mtry = 2, replace = TRUE, 
fraction = 0.632, minsplit = 120), data = BostonHousing, 
processors = 1, model = linearModel)

 rfout		## should produce descriptive output for mobForestOutput object
}
}
}
