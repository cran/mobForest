\name{treePredictions}
\alias{treePredictions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions from tree model
}
\description{
Predictions from tree model
}
\usage{
treePredictions(j, data, tree)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{j}{
jth observation
}
  \item{data}{
A data frame containing the variables in the model.
}
  \item{tree}{
An object of class mob inheriting from \code{\linkS4class{BinaryTree}}
}
}
\details{
This method computes predicted outcome for each observation in the data frame using the tree model supplied as an input argument.
}
\value{
A vector of predcited outcome
}