\name{mobForestOutput-class}
\Rdversion{1.1}
\docType{class}
\alias{mobForestOutput-class}
\alias{getPredictedValues,mobForestOutput,logical,logical-method}
\alias{getVarimp,mobForestOutput-method}
\alias{PredictiveAccuracy,mobForestOutput,logical,logical-method}
\alias{residualPlot,mobForestOutput-method}
\alias{show,mobForestOutput-method}
\alias{varimplot,mobForestOutput-method}

\title{Class \code{"mobForestOutput"}}
\description{
Random Forest output for model based recursive partitioning
}
\section{Objects from the Class}{
Objects can be created by \code{\link{mobForest_output}}.
}
\section{Slots}{
  \describe{
    \item{\code{oobPredictions}:}{Object of class \code{"predictionOutput"} for out-of-bag cases. }
    \item{\code{GeneralPredictions}:}{Object of class \code{"predictionOutput"} for all cases. }
    \item{\code{NewDataPredictions}:}{Object of class \code{"predictionOutput"} for new Test dataset cases.}
    \item{\code{VarimpObject}:}{Object of class \code{"varimpOutput"} ~~ }
    \item{\code{modelUsed}:}{ the model considered during random forest analysis. }
    \item{\code{fam}:}{ error distribution assumption }
    \item{\code{train.response}:}{data.frame of predictions for training set }
    \item{\code{new.response}:}{data.frame of predictions for new test dataset }
  }
}
\section{Methods}{
  \describe{
    \item{getPredictedValues}{\code{signature(object = "mobForestOutput", OOB = "logical", newdata = "logical")}: 
	get predictions summarized across trees for OOB cases or learning data or new test data. }
    \item{getVarimp}{\code{signature(object = "mobForestOutput")}: 
	get variable importance scores computed through random forest analysis }
    \item{PredictiveAccuracy}{\code{signature(object = "mobForestOutput", newdata = "logical", plot = "logical")}: 
	predictive performance across all trees	}
    \item{residualPlot}{\code{signature(object = "mobForestOutput")}: 
	produces two plots on same panel: a) histogram of residuals, b) predicted Vs residuals }
    \item{show}{\code{signature(object = "mobForestOutput")}: 
	print \code{mobForestOutput} object}
    \item{varimplot}{\code{signature(object = "mobForestOutput")}: 
	produces a plot with variable importance scores on X-axis and variable names on Y-axis.}
	 }
}
\seealso{
\code{\linkS4class{predictionOutput}},  \code{\linkS4class{varimpOutput}}
}
\examples{
\donttest{
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

## print method for mobForestOutput object
show(rfout)
}
}
}
\keyword{classes}
