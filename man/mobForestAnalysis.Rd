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
mobForestAnalysis(formula, partitionVariables, data, mobForest.controls = mobForest_control(), newTestData = as.data.frame(matrix(0, 0, 0)), processors = 1, model = linearModel, family = NULL, prob.cutoff = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
	An object of class formula specifying the model. This should be of type y ~ x1 + ... + xk, where the variables x1, x2, ..., xk are predictor variables and y represents an outcome variable. This model is referred to as the node model}
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
	A model of class \code{"\link[=StatModel-class]{StatModel}"} used for fitting observations in current node. This parameter allows fitting a linear model or generalized linear model with formula y ~ x1 + ... + xk. The Parameter "linearModel" fits linear model. The parameter "glinearModel" fits Poisson or logistic regression model depending upon the specification of parameter "family" (explained next). If "family"is specified as binomial() then logistic regression is performed. If the "family" is specified as poisson() then Poisson regression is performed.  }
  \item{family}{
	A description of error distribution and link function to be used in the model. This parameter needs to be specified if generalized linear model is considered. The parameter "binomial()" is to be specified when logistic regression is considered and "poisson()" when Poisson regression is considered as the node model. The values allowed for this parameter are binomial() and poisson().}
  \item{prob.cutoff}{
	In case of logistic regression as a node model, the predicted probabilities for OOB cases are converted into classes (yes/no, high/low, etc as specified) based on this probability cutoff. If logistic regression is not cosndiered as nodel model, the prob.cutoff = NULL. By default it is 0.5 when parameter not specieid (and logistic regression cponsidered). }
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
BostonHousing <- BostonHousing[1:90,c("rad", "tax", "crim", "medv", "lstat")]

## Recursive partitioning based on linear regression model medv ~ lstat with 3 trees. 
## 1 core/processor used. Supply more processors using 'processors' argument
rfout <- mobForestAnalysis(as.formula(medv ~ lstat), 
c("rad", "tax", "crim"), mobForest.controls = 
mobForest_control(ntree = 3, mtry = 2, replace = TRUE, 
alpha = 0.05, bonferroni = TRUE, minsplit = 25), 
data = BostonHousing, processors = 1, model = linearModel)

\donttest{ 
rfout		## should produce descriptive output for mobForestOutput object
}
  
}
}

