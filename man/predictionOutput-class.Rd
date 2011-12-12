\name{predictionOutput-class}
\Rdversion{1.1}
\docType{class}
\alias{predictionOutput-class}

\title{Class \code{"predictionOutput"}}
\description{
The object of this class stores predictions and predictive accuracy estimates.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{\link{prediction_output}}.

%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{predMat}:}{a matrix containing three columns: mean and standard deviation of predictions across trees and residuals (observed - predicted) }
    \item{\code{R2}:}{a vector of predictive accuracy estimates (ranging between 0 and 1) for each tree }
    \item{\code{mse}:}{a vector of mean-square-estimates (MSE) for each tree. Valid only if the outcome is continuous.}
    \item{\code{overallR2}:}{Overall predictive accuracy estimate obtained by combining predictions across trees.}
    \item{\code{predType}:}{character specifying the type of prediction output: out-of-bag cases, learning set, or new Test dataset}
  }
}
\section{Methods}{
No methods defined with class "predictionOutput" in the signature.
}

\seealso{
\code{\link{prediction_output}}, \code{\link{PredictiveAccuracy}}
}
\keyword{classes}
