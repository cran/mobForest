\name{logistic_accuracy}
\alias{logistic_accuracy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Contingency table: Predictive Vs Observed outcome
}
\description{
This function takes predicted probabilities (for out of bag cases) obtained through logistic regression-based tree models and converts them into binary classes (based on specified probability threshold). The predicted classifications are then compared to actual binary response.
}
\usage{
logistic_accuracy(response, predicted, prob.thresh)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{response}{
A vector of binary classes of out-of-cases for a given tree
}
  \item{predicted}{
A vector of predicted probabilities of out-of-cases using same tree
}
  \item{prob.thresh}{
Probibility threshold for classification. It is 0.5 by default.
}
}
