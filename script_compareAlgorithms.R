source("computeLogLikelihood.R")
source("siteAlgo.R")
source("masterAlgo.R")
source("centralizedNewtonRaphson.R")

# algorithm parameters
numIter = 100
tol = 10^(-10)

# initialization
bInitialization <-rep(0, times = dim(X)[2])
llInitialization = computeLogLikelihood(X,y,bInitialization)


# initialize recording arrays
bRecord = matrix(NA,nrow=numIter, ncol=dim(X)[2])
bCentralizedRecord = matrix(NA,nrow=numIter, ncol=dim(X)[2])
llRecord = matrix(NA,nrow=numIter, ncol=1)
llCentralizedRecord = matrix(NA,nrow=numIter, ncol=1)


bNew = bInitialization
llNew = llInitialization

bNewCentralized = bInitialization
llNewCentralized = llInitialization
cat("NR start\n")
for (i in 1:100)
{
  cat("Iteration ",i,"\n")
  # use new coefficients in next iteration
  b = bNew
  ll = llNew
  bCentralized = bNewCentralized
  llCentralized = llNewCentralized
  
  # compute NR components on each site
  output1 <- siteAlgo(X1,y1,b)
  output2 <- siteAlgo(X2,y2,b)
  output3 <- siteAlgo(X3,y3,b)
  output4 <- siteAlgo(X4,y4,b)
  
  # combine results from all sites
  firstTermCombined = output1[[1]] + output2[[1]] + output3[[1]] + output4[[1]]
  secondTermCombined = output1[[2]] + output2[[2]] + output3[[2]] + output4[[2]]
  
  # compute NR step
  bNew <- masterAlgo(firstTermCombined,secondTermCombined,b)
  # compute new log likelihood
  llNew <- computeLogLikelihood(X,y,bNew)
  
  # compute NR centralized (for comparison)
  bNewCentralized <- centralizedNewtonRaphson(X,y,bCentralized)
  llNewCentralized <- computeLogLikelihood(X,y,bNewCentralized)

  # record b and log likelihood
  bRecord[i,] = bNew
  bCentralizedRecord[i,] = bNewCentralized
  llRecord[i] = llNew
  llCentralizedRecord[i] = llNewCentralized
  if (sum(abs(b - bNew)) < tol) {break}
}
cat("NR end\n")

# fit coefficients using glm()
NormalModel   <- glm(HPV~., data=fit_data, family = binomial)
bGlm = NormalModel$coefficients
# compute log likelihood for glm()
llGlm <- computeLogLikelihood(X,y,bGlm)

# print log likelihoods to console
cat("Log likelihood\n")
cat("Distributed NR:",ll,"\n")
cat("Centralized NR:",llCentralized,"\n")
cat("glm():",llGlm,"\n")

# store coeffs in matrix for displaying in console
coeffTable = matrix(NA,nrow=3, ncol=dim(X)[2])
coeffTable[1,] = t(b)
coeffTable[2,] = t(bCentralized)
coeffTable[3,] = bGlm

# print coefficients to console
cat("Coefficients\n")
print(coeffTable)

# plot difference between distributed NR and glm()
bVsBGlmDiff = sweep(bRecord,2,bGlm)
plot(rowSums(abs(bVsBGlmDiff)))

