masterAlgo <-function(firstTermCombined,secondTermCombined,b)
{
  # This function aggregates the site-dependent terms for a single Newton-Raphson step for logistic regression. 
  # It assumes that the individual firstTerm and secondTerms from all sides have been summed up before passed into this function. 
  #Based on 
  # Jiang, W., Li, P., Wang, S., Wu, Y., Xue, M., Ohno-Machado, L. and Jiang, X., 2013. WebGLORE: a web service for Grid LOgistic REgression. Bioinformatics, 29(24), pp.3238-3240.
  
  firstTermInverse = solve(firstTermCombined)
  bNew = b + firstTermInverse %*% secondTermCombined
  
  return(bNew)
}