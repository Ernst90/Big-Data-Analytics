## Big Data Analytics (STATS 5016) - Masters Level 
## Project: Stochastic Optimization ALgorithms
## GUID: 2383746W


### Cost function and Gradient cost

## Cost function for linear regression
lmf <- function(y, X, theta, ndata) {
  sum((X%*%theta-y)^2)/(2*ndata)
}

## Gradient of cost function for linear regression
lmgrad <- function(y, X, theta, ndata){
  t(X)%*%(X%*%theta-y)/ndata
}
