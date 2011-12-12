\name{mobForest_output}
\alias{mobForest_output}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Model-based Random forest Object
}
\description{
Random Forest Output object that stores all the results including predictions, variable importance matrix, model, family of error distributions, and observed responses.
}
\usage{
mobForest_output(oobPredictions, GeneralPredictions, NewDataPredictions, 
VarimpObject, modelUsed, fam, train.response, 
new.response = data.frame(matrix(0, 0, 0)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{oobPredictions}{
Predictions on Out-of-bag data
}
  \item{GeneralPredictions}{
Predictions on learning data 
}
  \item{NewDataPredictions}{
Predictions on new test data 
}
  \item{VarimpObject}{
variable importance object
}
  \item{modelUsed}{
Model used
}
  \item{fam}{
a description of the error distribution and link function to be used in the model.
}
  \item{train.response}{
Response outcome of training data
}
  \item{new.response}{
Response outcome of test data
}
}

\seealso{
\code{\linkS4class{predictionOutput}}, \code{\linkS4class{varimpOutput}} 
}
