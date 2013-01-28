### mob Forest Analysis - all the classes are here
### this version includes mob_control class within mobForestControl class
### Version 3 deletes predictControl class
### Version 4 was made because parseFormula() had a bug that needed to be rectified. 
### Version 5 for adding the generalized linear model functionality

# setClass("mobTreeControl", representation(mob.control = "list"), prototype = list(alpha = 0.05, bonferroni = TRUE, minsplit = 20, trim = 0.1, objfun = deviance, breakties = FALSE, parm = NULL, verbose = FALSE))
# mobTree_control <- function(alpha = 0.05, bonferroni = TRUE, minsplit = 20, trim = 0.1, objfun = deviance, breakties = FALSE, parm = NULL, verbose = FALSE)
# {
#   mob.control = mob_control(alpha = alpha, bonferroni = bonferroni, minsplit = minsplit, trim = trim, objfun = objfun, breakties = breakties, parm = parm, verbose = verbose)
#   class(mob.control) <- "list"
#   rval <- new("mobTreeControl", mob.control = mob.control)
#   return(rval)
# }

setClass("mobForestControl", representation(ntree = "numeric", mtry = "numeric", replace = "logical", fraction="numeric", mob.control = "list"), prototype = list(ntree = 300, mtry = 0, replace = FALSE, fraction = 0.632))
mobForest_control <- function(ntree = 300, mtry = 0, replace = FALSE, fraction = 0.632, alpha = 1, bonferroni = FALSE, minsplit = 20, trim = 0.1, objfun = deviance, breakties = FALSE, parm = NULL, verbose = FALSE)
{
  mob.control = mob_control(alpha = alpha, bonferroni = bonferroni, minsplit = minsplit, trim = trim, objfun = objfun, breakties = breakties, parm = parm, verbose = verbose)
  class(mob.control) <- "list"
  rval <- new("mobForestControl", ntree = ntree, mtry = mtry, replace = replace, fraction = fraction, mob.control = mob.control)
  return(rval)
}

setClass("predictionOutput", representation(predMat = "matrix", R2 = "numeric", mse = "numeric", overallR2 = "numeric", predType = "character"), prototype = list(predMat = matrix(0,0,0), R2 = numeric(), mse = numeric(), overallR2 = numeric(), predType = character()))
prediction_output <- function(predMean = numeric(), predSd = numeric(), residual = numeric(), R2 = numeric(), mse = numeric(), overallR2 = numeric(), predType = character())
{	
	predMat = matrix(0, nrow=length(predMean), ncol=3)
	predMat[,1] = predMean
	predMat[,2] = predSd
	predMat[,3] = residual
	colnames(predMat) = c("PredMean", "PredStdev", "Residual")
	rval <- new("predictionOutput", predMat = predMat, R2 = R2, mse = mse, overallR2 = overallR2, predType = predType)
	return(rval)
}

setClass("varimpOutput", representation(varimpMatrix = "matrix"), prototype = list(varimpMatrix = matrix(0,0,0)))
varimp_output <- function(varimpMatrix)
{	
	rval <- new("varimpOutput", varimpMatrix = varimpMatrix)    
    return(rval)
}

setClass("mobForestOutput", representation(oobPredictions = "predictionOutput", GeneralPredictions = "predictionOutput", NewDataPredictions = "predictionOutput", VarimpObject = "varimpOutput", modelUsed = "character", fam = "character", train.response = "data.frame", new.response = "data.frame"))
mobForest_output <- function(oobPredictions, GeneralPredictions, NewDataPredictions, VarimpObject, modelUsed, fam, train.response, new.response = data.frame(matrix(0,0,0)))
{
	rval <- new("mobForestOutput",oobPredictions = oobPredictions, GeneralPredictions = GeneralPredictions, NewDataPredictions = NewDataPredictions, VarimpObject = VarimpObject, modelUsed = modelUsed, fam = fam, train.response = train.response, new.response = new.response)
    return(rval)
}

setMethod("show", "mobForestOutput", function(object) {
		
	rf <- object
	cat("\tRandom Forest of Model Based Recursive Partitioning Trees\n\n")
	cat(paste("Number of trees:",ncol((rf@VarimpObject)@varimpMatrix),"\n\n"))
	cat(paste("Model Used:",rf@modelUsed,"\n\n"))	 
})

setGeneric("getVarimp", function(object) standardGeneric("getVarimp"))
setMethod("getVarimp", signature(object="mobForestOutput"), function(object) {
	
	rf <- object
	varImp.scores <- apply((rf@VarimpObject)@varimpMatrix, 1, mean, na.rm=T)	
	return(sort(varImp.scores, decreasing=T))
})

setGeneric("varimplot", function(object) standardGeneric("varimplot"))
setMethod("varimplot", signature(object="mobForestOutput"), function(object) {
	
	rf <- object
	library(lattice)
	varImp.scores <- apply((rf@VarimpObject)@varimpMatrix, 1, mean, na.rm=T)	
	par(mfrow = c(1,2))
	#boxplot(t(rf@VarimpObject@varimpMatrix))
	lattice::dotplot(sort(varImp.scores),xlab="Variable Importance in the data",
        panel = function(x,y){
           panel.dotplot(x,y,col='darkblue', pch=16, cex=1.1, main="Variance Importance Plot")
           panel.abline(v=abs(min(varImp.scores)), col='red',
           lty='longdash', lwd=2)
           panel.abline(v=0, col='blue')
})	
})

setGeneric("getPredictedValues", function(object, OOB = TRUE, newdata = FALSE) standardGeneric("getPredictedValues"))
setMethod("getPredictedValues", signature(object="mobForestOutput", OOB="ANY", newdata="ANY"), function(object, OOB, newdata) {
	
	#if(missing(object)) cat("This function expects object of 'mobForestOutput', returned by mobForestAnalysis(), as its first argument\n")
	rf <- object	
	if (nrow(rf@NewDataPredictions@predMat) == 0 && newdata == TRUE)
		stop("Predicted values were only computed on original data. Please set newdata = FALSE and run the getPredictedValues() again. Or you can re-run the mobForestAnalysis() with 'newTestdata' parameter not missing and later use getPredictedValues() to get predicted values on the new test data.")
	rval <- c()
	if(newdata == FALSE)
	{
		if(OOB == TRUE)
		{
			rval <- (rf@oobPredictions)@predMat			
		} else {
			rval <- (rf@GeneralPredictions)@predMat
		}
	} else {
		rval <- (rf@NewDataPredictions)@predMat		
	}		
	return(rval)
})

setGeneric("residualPlot", function(object) standardGeneric("residualPlot"))
setMethod("residualPlot", signature(object="mobForestOutput"), function(object) {
  	
	rf <- object
  if(rf@fam == "binomial" | rf@fam == "poisson")
  {
    cat("Residual Plot not produced when logistic of Poisson regression is considered as the node model\n")
    #return;
  } else {
		par(mfrow = c(2,1))
	  plot(rf@oobPredictions@predMat[,1], rf@oobPredictions@predMat[,3], xlab="Out-of-bag predcitions", ylab="Out-of-bag residuals")
	  hist(rf@oobPredictions@predMat[,3], main="Out-of-bag residuals histogram", xlab="Out-of-bag residuals", breaks=50)
  }  
})

setGeneric("PredictiveAccuracy", function(object, newdata = FALSE, prob.cutoff = NULL, plot = TRUE) standardGeneric("PredictiveAccuracy"))
setMethod("PredictiveAccuracy", signature(object="mobForestOutput", newdata="ANY", prob.cutoff="ANY", plot="ANY"), function(object, newdata, prob.cutoff, plot) {
	
	rf <- object
	rval <- list()
	vecR2 = NULL
	vecMSE = NULL
	if(newdata == FALSE)
	{		
		ss = nrow((rf@GeneralPredictions)@predMat)
		rval <- list((rf@oobPredictions)@R2, (rf@oobPredictions)@mse, (rf@oobPredictions)@overallR2, sum((rf@oobPredictions)@predMat[,3]**2, na.rm=T)/ss, (rf@GeneralPredictions)@R2, (rf@GeneralPredictions)@mse, (rf@GeneralPredictions)@overallR2, sum((rf@GeneralPredictions)@predMat[,3]**2, na.rm=T)/ss, rf@modelUsed, rf@fam, prob.cutoff)
		names(rval) <- c("oob.R2", "oob.mse", "oob.OverallR2", "oob.OverallMSE", "General.R2", "General.mse", "General.OverallR2", "General.OverallMSE", "modelUsed", "fam", "prob.cutoff")
		vecR2 = c(rval$oob.OverallR2, min(rval$oob.R2), max(rval$oob.R2))
		vecMSE = c(rval$oob.OverallMSE, min(rval$oob.mse), max(rval$oob.mse))
	} else {		
		ss1 = nrow((rf@GeneralPredictions)@predMat)
		ss2 = nrow((rf@NewDataPredictions)@predMat)
		rval <- list((rf@oobPredictions)@R2, (rf@oobPredictions)@mse, (rf@oobPredictions)@overallR2, sum((rf@oobPredictions)@predMat[,3]**2, na.rm=T)/ss1, (rf@GeneralPredictions)@R2, (rf@GeneralPredictions)@mse, (rf@GeneralPredictions)@overallR2, sum((rf@GeneralPredictions)@predMat[,3]**2, na.rm=T)/ss1, rf@modelUsed, rf@fam, (rf@NewDataPredictions)@R2, (rf@NewDataPredictions)@overallR2, sum((rf@NewDataPredictions)@predMat[,3]**2, na.rm=T)/ss2, prob.cutoff)
		names(rval) <- c("oob.R2", "oob.mse", "oob.OverallR2", "oob.OverallMSE", "General.R2", "General.mse", "General.OverallR2", "General.OverallMSE", "modelUsed", "fam", "Newdata.R2", "Newdata.OverallR2", "Newdata.OverallMSE", "prob.cutoff")
		vecR2 = c(rval$oob.OverallR2, min(rval$oob.R2), max(rval$oob.R2), rval$Newdata.OverallR2)
		vecMSE = c(rval$oob.OverallMSE, min(rval$oob.mse), max(rval$oob.mse), rval$Newdata.OverallMSE)
	}
	
	newd = which(regexpr("Newdata", names(rval)) > 0)
	if(rf@fam == "binomial")
	{
		rval$train.response = rf@train.response
		rval$new.response = rf@new.response
		rval$oobPredMean = rf@oobPredictions@predMat[,1]
		rval$genPredMean = rf@GeneralPredictions@predMat[,1]
		#rval <- list(rval, rf@train.response, rf@new.response, rf@oobPredictions@predMat[,1], rf@GeneralPredictions@predMat[,1])
		#names(rval) <- c(names(rval), "train.response", "new.response", "oobPredMean", "genPredMean")
		#if(length(newd) > 0) { rval <- list(rval, rf@NewDataPredictions@predMat[,1]); names(rval) <- c(names(rval), "newPredMean") }
		if(length(newd) > 0) rval$newPredMean <- rf@NewDataPredictions@predMat[,1]
	}
	
	if(plot == TRUE)
	{		
		xlab = ""
		r = 2		
		if(rf@fam == "binomial")
		{
			r = r - 1
			xlab = "Proportion of subjects correctly classified"			
		} else {
			xlab = expression(R**2)
		}		
		
		xlim1 = ifelse(min(vecR2) < 0, 0, min(vecR2))
		#xlim2 = max(vecR2) + (max(vecR2) - xlim1)/5		
		xlim2 = ifelse(max(vecR2) > 1, 1, max(vecR2))
		if(length(newd) == 0) par(mfrow = c(r,2))
		if(length(newd) > 0) par(mfrow = c((r + 1),2))		
		hist(rval$oob.R2, main="OOB performance (Tree Level)", xlab=xlab, xlim=c(xlim1, xlim2))
		box()
		plot(c(rval$oob.OverallR2, rval$oob.OverallR2), c(1,2), axes=F, type="l", col="red", lwd=2, main = "OOB performance (Forest Level)", xlab=xlab, ylab="", xlim=c(xlim1, xlim2))
		axis(1)		
		box()
		if(rf@fam != "binomial") 
		{
			hist(rval$oob.mse, xlab="MSE", main = "OOB performance (Tree Level)", xlim=c(min(vecMSE), max(vecMSE)))
			box()
			plot(c(rval$oob.OverallMSE, rval$oob.OverallMSE), c(1,2), axes=F, type="l", col="red", lwd=2, xlab="MSE", ylab="", main = "OOB performance (Forest Level)", xlim=c(min(vecMSE), max(vecMSE)))
			axis(1)			
			box()
		}
		if(length(newd) > 0)
		{
			if(length(rval$Newdata.OverallR2) != 0) 
			{
				plot(c(rval$Newdata.OverallR2,rval$Newdata.OverallR2), c(1,2), axes=F, type="l", col="red", lwd=2, xlab=xlab, ylab="", xlim=c(xlim1, xlim2), main="Validation Performance")
				axis(1)
				box()
				if(rf@fam != "binomial")
				{
					plot(c(rval$Newdata.OverallMSE,rval$Newdata.OverallMSE), c(1,2), axes=F, type="l", col="red", lwd=2, xlab="MSE", ylab="", main="Validation Performance", xlim=c(min(vecMSE), max(vecMSE)))
					axis(1)
					box()
				}
			} else {
				cat("No Performance Plot for new test data because 'newTestData' argument was not supplied to mobForestAnalysis() function \n")
			}
		}
		#title("Predictive Accuracy", outer = TRUE, line = -1)
	}
	class(rval) <- "predAccuracyEstimates"	
	return(rval)
})

print.predAccuracyEstimates <- function(x,...)
{
	pacc <- x
	pred.accuracy <- x
	cat("MobForest Predictive Accuracy Report:\n\n")
	cat("Tree Level Summary \n\n")
	if(pred.accuracy$fam != "binomial") cat("Pseudo R2 Report \n")
	if(pred.accuracy$fam == "binomial") cat("Proportion of subjects correctly classified\n")
	cat("Data Used\tMin.\tQ1\tMedian\tMean\tQ3\tMax.\n")
	oobR2.qs = as.numeric(quantile(pacc$oob.R2, probs = c(0, 0.25, 0.5, 0.75, 1)))
	cat("OOB (Training) ", round(oobR2.qs[1:3],3), round(mean(pacc$oob.R2),3), round(oobR2.qs[4:5],3), sep="\t")
	cat("\n")
	genR2.qs = as.numeric(quantile(pacc$General.R2, probs = c(0, 0.25, 0.5, 0.75, 1)))
	cat("All (Training)", round(genR2.qs[1:3],3), round(mean(pacc$General.R2),3), round(genR2.qs[4:5],3), sep="\t")
	cat("\n\n")
	
	if(pred.accuracy$fam != "binomial")
	{
		cat("MSE Report \n")
		cat("Data Used\tMin.\tQ1\tMedian\tMean\tQ3\tMax.\n")
		oobMSE.qs = as.numeric(quantile(pacc$oob.mse, probs = c(0, 0.25, 0.5, 0.75, 1)))
		cat("OOB (Training)", round(oobMSE.qs[1:3],3), round(mean(pacc$oob.mse),3), round(oobMSE.qs[4:5],3), sep="\t")
		cat("\n")
		genMSE.qs = as.numeric(quantile(pacc$General.mse, probs = c(0, 0.25, 0.5, 0.75, 1)))
		cat("All (Training)", round(genMSE.qs[1:3],3), round(mean(pacc$General.mse),3), round(genMSE.qs[4:5],3), sep="\t")
		cat("\n\n")
	}
	
	cat("Forest Level Summary \n\n")	
	if(pred.accuracy$fam == "binomial")
	{		
		cat("OOB (Training) Data: \n\n")
		logistic_accuracy(pacc$train.response, pacc$oobPredMean, pacc$prob.cutoff)	
		cat("All (Training) Data: \n\n")
		logistic_accuracy(pacc$train.response, pacc$genPredMean, pacc$prob.cutoff)
		if(!is.null(pacc$Newdata.OverallR2))
		{
			cat("Validation Data: \n\n")
			logistic_accuracy(pacc$new.response, pacc$newPredMean, pacc$prob.cutoff)
		}		
	} else {
		cat("Data Used\tPseudo R2\tMSE\n") 
		cat("OOB (Training)", formatC(c(pacc$oob.OverallR2, pacc$oob.OverallMSE), format="f", digits=7, drop0trailing=FALSE), sep="\t")
		cat("\n")
		cat("All (Training)", formatC(c(pacc$General.OverallR2, pacc$General.OverallMSE), format="f", digits=7, drop0trailing=FALSE), sep="\t")
		cat("\n")
		if(!is.null(pacc$Newdata.OverallR2))
		{
			if (pred.accuracy$fam != "binomial") cat("Validation", formatC(c(pacc$Newdata.OverallR2, pacc$Newdata.OverallMSE), format="f", digits=7, drop0trailing=FALSE), sep="\t")
			cat("\n")		
		}
	}
}

logistic_accuracy <- function(response, predicted, prob.thresh)
{
	if(is.null(prob.thresh)) prob.thresh = 0.5
  PredClass = rep(levels(response[,1])[1], length(predicted))
	PredClass[which(predicted > prob.thresh)] = levels(response[,1])[2]
	Tab = table(Data_Class = response[,1], Predicted_Class = PredClass)
	#rownames(Tab) = c(paste("True", levels(response[,1])[1]), paste("True", levels(response[,1])[2]))
	#colnames(Tab) = c(paste("Predicted", levels(predicted)[1]), paste("Predicted", levels(predicted)[2]))
	print(Tab)
	cat("\n\n")
}
