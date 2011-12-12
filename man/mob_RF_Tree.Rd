\name{mob_RF_Tree}
\alias{mob_RF_Tree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
model based recursive partitioning - randomized subset of partition variables considered during each split.
}
\description{
model based recursive partitioning - randomized subset of partition variables considered during each split.
}
\usage{
mob_RF_Tree(mainModel, partitionVars, mtry, weights, data = list(), 
na.action = na.omit, model = glinearModel, control = mob_control(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mainModel}{
A model in character format
}
  \item{partitionVars}{
A vector of partition variables
}
  \item{mtry}{
A Random subset of partition variables to be considered at each node of decesion tree
}
  \item{weights}{
An optional vector of weights, as described in \link[=mob]{mob}
}
  \item{data}{
A data frame containing the variables in the model.
}
  \item{na.action}{
A function which indicates what should happen when the data contain NAs, as described in \link[=mob]{mob}
}
  \item{model}{
A model of class \code{\linkS4class{StatModel}}
}
  \item{control}{
A list with control parameters as returned by \link[=mob_control]{mob_control}
}
  \item{\dots}{
Additional arguments passed to the fit call for the model.
}
}
\details{
The \link[=mob]{mob} function in party package is modified so that a random subset of predictor variables are considered during each split. mtry represents the number of predictor variables to be considered during each split. 
}
\value{
An object of class mob inheriting from \code{\linkS4class{BinaryTree}}. Every node of the tree is additionally associated with a fitted model.
}
\references{
Achim Zeileis, Torsten Hothorn, and Kurt Hornik (2008). Model-Based Recursive Partitioning. \emph{Journal of Computational and Graphical Statistics}, 17(2), 492-514.\cr
}