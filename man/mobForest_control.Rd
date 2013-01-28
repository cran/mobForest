\name{mobForest_control}
\alias{mobForest_control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Control parameters for random forest 
}
\description{
Various parameters that control the forest growing.
}
\usage{
mobForest_control(ntree = 300, mtry = 0, replace = FALSE, fraction = 0.632, 
alpha = 1, bonferroni = FALSE, minsplit = 20, trim = 0.1, objfun = deviance, 
breakties = FALSE, parm = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ntree}{
Number of trees to be constructed in forest (default = 300)
}
  \item{mtry}{
number of input variables randomly sampled as candidates at each node
}
  \item{replace}{
TRUE or FALSE. replace = TRUE (default) performs bootstrapping. replace = FALSE performs sampling without replacement.
}
  \item{fraction}{
number of observations to draw without replacement (only relevant if replace = FALSE)
}
  \item{alpha}{
	A node is considered for splitting if the p value for any partitioning variable in that node falls below alpha (default 0.05) as mentioned in \link[=mob_control]{mob_control()}
}
  \item{bonferroni}{
	logical, Should p values be Bonferroni corrected? (default TRUE) as mentioned in \link[=mob_control]{mob_control()}
}
  \item{minsplit}{
	integer, The minimum number of observations in a node (default 20) as mentioned in \link[=mob_control]{mob_control()}
}
  \item{trim}{
	numeric, as defined in \link[=mob_control]{mob_control()}
}
  \item{objfun}{
	function, as defined in \link[=mob_control]{mob_control()}
}
  \item{breakties}{
	logical, as defined in \link[=mob_control]{mob_control()}
}
  \item{parm}{
	numeric or vector, as defined in \link[=mob_control]{mob_control()}
}
  \item{verbose}{
	logical, as defined in \link[=mob_control]{mob_control()}
}
}
\details{
This function is used to set up forest controls. The mob_control (from party 'package') object is used to set up control parameters for single tree model.
}
\value{
An object of class \code{\linkS4class{mobForestControl}}.
}
\references{
Achim Zeileis, Torsten Hothorn, and Kurt Hornik (2008). Model-Based Recursive Partitioning. \emph{Journal of Computational and Graphical Statistics}, 17(2), 492-514.\cr
}
\examples{
## create forest controls before starting random forest analysis
mobForest.controls = mobForest_control(ntree = 400, mtry = 4, 
replace = TRUE, minsplit = 200)
}
