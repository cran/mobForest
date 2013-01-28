\name{PredictiveAccuracy}
\alias{PredictiveAccuracy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictive Performance
}
\description{
Predictive performance across all trees
}
\usage{
PredictiveAccuracy(object, newdata = FALSE, prob.cutoff = NULL, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{mobForestOutput}
}
  \item{newdata}{
A logical value specifying if the performance needs to be summarized on test data supplied as \code{newTestdata} argument to \code{mobForestAnalysis} function.
}
  \item{prob.cutoff}{
Predicted probabilities converted into classes (Yes/No, 1/0) based on this probability threshold. Only used for producing predicted Vs actual classes table.
}
  \item{plot}{
A logical value specifying if the use wishes to veiw performance plots
}
}
\value{
A list with performance parameters with class "predAccuracyEstimates"
%%  \item{oob.R2 }{A vector of predictive accuracy estimates (ranging between 0 and 1) measured on Out-of-bag cases for each tree}
%%  \item{oob.mse }{A vector of MSE for Out-of-bag data for each tree. Valid only if the outcome is continuous.}
%%  \item{oob.OverallR2 }{ Overall predictive accuracy measured by combining Out-of-bag predictions across trees.}
%%  \item{oob.OverallMSE }{ Overall MSE measured by combining Out-of-bag predictions across trees.}
%%  \item{General.R2 }{ A vector of predictive accuracies (ranging between 0 and 1) measured on complete learning data for each tree }
%%  \item{General.mse }{ A vector of MSE measured on complete learning data for each tree. Valid only if the outcome is continuous. }
%%  \item{General.OverallR2 }{ Overall predictive accuracy measured by combining predictions across trees. }
%%  \item{General.OverallMSE }{ Overall MSE measured by combining predictions across trees. Valid only if the outcome is continuous.}
%%  \item{modelUsed }{ The node model and partition valriables used for analysis}
%%  \item{family }{ Error distribution assumptions of the model}
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

## get predictive performance estimates and produce a performance plot 
pacc <- PredictiveAccuracy(rfout)
}
