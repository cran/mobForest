\name{mobForestAnalysis}
\alias{mobForestAnalysis}
\alias{bootstrap}
\alias{getmobForestObject.GLM}
\alias{getmobForestObject.LM}
\title{ Model-based Random Forest Analysis }
\description{
 Main function that takes all the neccessary arguments to start model-based random forest analysis.
}
\usage{
mobForestAnalysis(formula, partitionVariables, data, mobForest.controls = 
mobForest_control(), newTestData = as.data.frame(matrix(0, 0, 0)), processors = 1, 
model = linearModel, family = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
	An object of class formula specifying the model. This should be of type y ~ x1 + ... + xk, where the variables x1, x2, ..., xk are predictor variables and y represents an outcome variable.}
  \item{partitionVariables}{
	A character vector specifying the partition variables }
  \item{data}{
	An input dataset that is used for constructing trees in random forest.}
  \item{mobForest.controls}{
	An object of class \code{"\linkS4class{mobForestControl}"} returned by \link[=mobForest_control]{mobForest_control()}, that contains parameters controlling the construction of random forest. }
  \item{newTestData}{
	A data frame representing test data for validating random forest model. This data is not used in in tree building process. }
  \item{processors}{
	A number of processors/cores on your computer that should be used for parallel computation. }
  \item{model}{
	A model of class \code{"\link[=StatModel-class]{StatModel}"} used for fitting observations in current node. This parameter allows fitting a linear model or generalized linear model with formula y ~ x1 + ... + xk. }
  \item{family}{
	A description of error distribution and link function to be used in the model. This parameter needs to be specified if generalized linear model is considered. }
}

\details{ 
\code{mobForestAnalysis} is the main function that takes all the input parameters - model, partition variables, and forest control parameters - and starts the model-based random forest analysis. \code{mobForestAnalysis} calls \code{bootstrap} function which constructs decesion trees, computes out-of-bag (OOB) predictions, OOB predictive 	accuracy and perturbation in OOB predictive accuracy through permutation. \code{bootstrap} constructs trees on multiple cores/processors simultaneously through parallel computation. Later, the \code{getmobForestObject} function wraps the analysis output into \code{mobForestOutput} object.\cr

Predictive accuracy estimates are computed using pseudo-R2 metric, defined as the proportion of total variation in outcome variable explained by a tree model on out-of-bag cases. R2 ranges from 0 to 1. R2 of zero suggests worst tree model (in terms of predicting outcome) and R2 of 1 suggests perfect tree model.\cr 
}
\value{
An object of class \code{\linkS4class{mobForestOutput}}.
}
\references{
Achim Zeileis, Torsten Hothorn, and Kurt Hornik (2008). Model-Based Recursive Partitioning. \emph{Journal of Computational and Graphical Statistics}, 17(2), 492-514.\cr

Hothorn, T., Hornik, K. and Zeileis, A. (2006) Unbiased recursive partitioning: A conditional inference framework, \emph{J Comput Graph Stat}, 15, 651-674.\cr

Strobl, C., Malley, J. and Tutz, G. (2009) An introduction to recursive partitioning: rationale, application, and characteristics of classification and regression trees, bagging, and random forests, \emph{Psychol Methods}, 14, 323-348.\cr
}

\seealso{
\link[=mobForest_control]{mobForest_control()}, \code{\link{mobForestOutput-class}}
}
\examples{
set.seed(290875)

if(require("mlbench")) {

## Random Forest analysis of model based recursive partitioning
## load data  
data("BostonHousing", package = "mlbench")

## Recursive partitioning based on linear regression model medv ~ lstat + rm with 8 trees. 
## 1 core/processor used. Supply more processors using 'processors' argument
rfout <- mobForestAnalysis(as.formula(medv ~ lstat + rm), c("zn", "indus", "chas", "nox", 
"age", "dis", "rad", "tax", "crim", "b", "ptratio"), mobForest.controls = 
mobForest_control(ntree = 8, mtry = 4, replace = TRUE, fraction = 0.632, minsplit = 80),
data = BostonHousing, 
processors = 1,
model = linearModel)

rfout		## should produce descriptive output for mobForestOutput object
  
}
}

