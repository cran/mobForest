\name{mobForestControl-class}
\Rdversion{1.1}
\docType{class}
\alias{mobForestControl-class}

\title{Class \code{"mobForestControl"}}
\description{
Control parameters for random forest
}
\section{Objects from the Class}{
Objects can be created by \code{\link{mobForest_control}}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{ntree}:}{number fo trees in forest}
    \item{\code{mtry}:}{random subset of partition variables to be considered at each node of decesion tree}
    \item{\code{replace}:}{ sampling with replacement or without replacement}
    \item{\code{fraction}:}{proportion of  observations to be sampled without replacement}
    \item{\code{mob.control}:}{An object of class \code{\link{mob_control}}, which can be obtained using \code{\link{mobForest_control}}}
  }
}
\section{Methods}{
No methods defined with class "mobForestControl" in the signature.
}
\examples{
#showClass("mobForestControl")
## The following code creates following forest controls: 400 trees to be 
## constructed, sampling with replacement, a node contains at least 200 observations
mobForest.controls = mobForest_control(ntree = 400, mtry = 4, replace = TRUE, 
minsplit = 200)
}
\keyword{classes}
