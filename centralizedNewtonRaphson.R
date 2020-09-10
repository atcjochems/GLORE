centralizedNewtonRaphson <-function(X,y,b)
{
  # This function computes a single Netwon-Raphson step for logistic regression.
  # Based on derivations presented by 
  # Jiang, W., Li, P., Wang, S., Wu, Y., Xue, M., Ohno-Machado, L. and Jiang, X., 2013. WebGLORE: a web service for Grid LOgistic REgression. Bioinformatics, 29(24), pp.3238-3240.
  
  
  # it is assumed that the feature matrix X contains a row of ones in the last column, which needs to be moved to the front
  constantAndFeatures = (cbind(X$c,X[,1:dim(X)[2]-1]))
  
  constantAndFeatures = sapply(constantAndFeatures,as.numeric)
  predictions = 1/(1 + exp(-(constantAndFeatures%*%as.numeric(t(b)))))
  predictions = as.numeric(predictions)
  
  # calculate first term 
  probabilityProduct = predictions * (1-predictions)
  W = diag(probabilityProduct)
  firstTerm = t(constantAndFeatures) %*% W %*% constantAndFeatures
  
  
  # calculate second term
  predictionError = y-predictions
  secondTerm = t(constantAndFeatures) %*% predictionError
  
  firstTermInverse = solve(firstTerm)
  bNew = b + firstTermInverse %*% secondTerm
  return(bNew)
  
}