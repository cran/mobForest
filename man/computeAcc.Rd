\name{computeAcc}
\alias{computeAcc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictive accuracy estimates across trees for logistic regression model
}
\description{
Computes predictive accuracy of response variable with binary outcome. The function takes observed and predicted binary responses as input arguments and computes proportion of observations classified in same group.
}
\usage{
computeAcc(response, predictions)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{
a vector of binary outcome 
}
  \item{predictions}{
a matrix of predicted probabilities (logit model) for out-of-bag observations for each tree
}
}

\value{
{ Predictive accuracy estimates ranging between 0 and 1 for each tree. Zero represents worst accuracy and 1 suggests best accuracy.}
}
