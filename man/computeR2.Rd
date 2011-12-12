\name{computeR2}
\alias{computeR2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictive accuracy estimates across trees for linear or poisson regression model 
}
\description{
pseudo R-square (R2) computation - proportion of total variance in response variable explained by the tree model. The function takes observed and predicted responses as input arguments and computes pseudo-R2 to determine how well the tree model fits the given data.
}
\usage{
computeR2(response, predictions)
}

\arguments{
  \item{response}{
a vector of actual response of outcome variable.
}
  \item{predictions}{
a vector of predicted response for the same outcome variable.
}
}

\value{
{ Predictive accuracy estimates ranging between 0 and 1. Zero represents worst accuracy and 1 suggests best accuracy.}
}
