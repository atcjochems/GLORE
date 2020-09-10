computeLogLikelihood <- function(X,y,coeffs)
{constantAndFeatures = cbind(X$c,X[,1:dim(X)[2]-1])

predictions <- 1/(1 + exp(-(as.matrix(constantAndFeatures)%*%as.numeric(t(coeffs)))))
logLikelihood <- sum(y*log(predictions) + (1-y)*log(1-predictions))
return(logLikelihood)
}
