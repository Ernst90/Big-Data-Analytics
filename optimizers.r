## Big Data Analytics (STATS 5016) - Masters Level 
## Project: Stochastic Optimization
## GUID: 2383746W


### Algorithm 1: Gradient descent (GD) given a cost function f and its gradient 
gd <- function(f, grad, y, X, theta0, npars, ndata, a, niters) {
  theta <- matrix(data=NA, nrow=niters, ncol=npars)
  cost <- vector(mode="numeric", length=niters)
  
  # Starting values
  theta[1, ] <- theta0
  cost[1] <- f(y, X, theta0, ndata)
  
  # Updating steps 
  for (i in 2:niters) {
    theta[i, ] <- theta[i-1, ]-a*grad(y, X, theta[i-1, ], ndata)
    cost[i] <- f(y, X, theta[i, ], ndata)
  }
  return(list(theta=theta, cost=cost))
}


### Algorithm 2: Stochastic gradient descent (SGD) given a cost function f and its gradient 
sgd <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples) {
  theta <- matrix(data=NA, nrow=niters, ncol=npars)
  cost <- vector(mode="numeric", length=niters)
  
  # Starting values
  theta[1, ] <- theta0
  cost[1] <- f(y, X, theta0, ndata)
  
  # Updating steps 
  for (i in 2:niters) {
    dataset <- cbind(y, X)
    random.index <- dataset[sample(nrow(dataset), nsubsamples), ]
    X.samples <- cbind(rep(1, times=nsubsamples), random.index[ ,2:npars+1])
    y.samples <- random.index[ ,1]
    theta[i, ] <- theta[i-1, ] - a*grad(y.samples, X.samples, theta[i-1, ], nsubsamples)
    cost[i] <- f(y.samples, X.samples, theta[i, ], nsubsamples)
  }
  return(list(theta=theta, cost=cost))
}


### Algorithm 3: Stochastic gradient descent with momentum (MSGD) given a cost function f and its gradient 
msgd <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, m0) {
  theta <- matrix(data=NA, nrow=niters, ncol=npars)
  m <- vector(mode="numeric", length=niters)
  cost <- vector(mode="numeric", length=niters)
  
  # Starting values
  theta[1, ] <- theta0
  m <- m0
  cost[1] <- f(y, X, theta0, ndata)
  
  # Updating steps 
  for (i in 2:niters) {
    dataset <- cbind(y, X)
    random.index <- dataset[sample(nrow(dataset), nsubsamples), ]
    X.samples <- cbind(rep(1, times=nsubsamples), random.index[ ,2:npars+1])
    y.samples <- random.index[ ,1]
    m <- b*m + (1-b)*grad(y.samples, X.samples, theta[i-1, ], nsubsamples)
    theta[i, ] <- theta[i-1, ]- a*m
    cost[i] <- f(y.samples, X.samples, theta[i, ], nsubsamples)
  }
  return(list(theta=theta, cost=cost))
}


### Algorithm 4: Stochastic gradient descent with Nesterov accelerated gradient (NAGSGD) given a cost function f and its gradient 
nagsgd <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, m0) {
  theta <- matrix(data=NA, nrow=niters, ncol=npars)
  m <- vector(mode="numeric", length=niters)
  cost <- vector(mode="numeric", length=niters)
  
  # Starting values
  theta[1, ] <- theta0
  m <- m0
  cost[1] <- f(y, X, theta0, ndata)
  
  # Updating steps 
  for (i in 2:niters) {
    dataset <- cbind(y, X)
    random.index <- dataset[sample(nrow(dataset), nsubsamples), ]
    X.samples <- cbind(rep(1, times=nsubsamples), random.index[ ,2:npars+1])
    y.samples <- random.index[ ,1]
    m <- b*m + (1-b)*grad(y.samples, X.samples, theta[i-1, ]- a*b*m, nsubsamples)
    theta[i, ] <- theta[i-1, ]- a*m
    cost[i] <- f(y.samples, X.samples, theta[i, ], nsubsamples)
  }
  return(list(theta=theta, cost=cost))
}


### Algorithm 5: AdaGrad given a cost function f and its gradient 
adagrad <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples, epsilon, G0) {
  theta <- matrix(data=NA, nrow=niters, ncol=npars)
  cost <- vector(mode="numeric", length=niters)
  G <- vector(mode="numeric", length=niters)
  
  # Starting values
  theta[1, ]<- theta0
  cost[1] <- f(y, X, theta0, ndata)
  G <- G0
  
  # Updating steps 
  for(i in 2:niters){
    dataset <- cbind(y, X)
    random.index <- dataset[sample(nrow(dataset), nsubsamples), ]
    X.samples <- cbind(rep(1, times=nsubsamples), random.index[ ,2:npars+1])
    y.samples <- random.index[ ,1]
    g <- grad(y.samples, X.samples, theta[i-1, ], nsubsamples)
    G <- G + g^2
    theta[i, ] <- theta[i-1, ] - a*g / (sqrt(G + epsilon))
    cost[i] <- f(y.samples, X.samples, theta[i, ], nsubsamples)
  }
  return(list(theta=theta, cost=cost))
}


### Algorithm 6: RMSProp given a cost function f and its gradient 
rmsprop <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples, c, epsilon, v0) {
  theta <- matrix(data=NA, nrow=niters, ncol=npars)
  cost <- vector(mode="numeric", length=niters)
  v <- vector(mode="numeric", length=niters)

  # Starting values
  theta[1, ]<- theta0
  cost[1] <- f(y, X, theta0, ndata)
  v <- v0 
  
  # Updating steps 
  for(i in 2:niters){
    dataset <- cbind(y, X)
    random.index <- dataset[sample(nrow(dataset), nsubsamples), ]
    X.samples <- cbind(rep(1, times=nsubsamples), random.index[ ,2:npars+1])
    y.samples <- random.index[ ,1]
    g <- grad(y.samples, X.samples, theta[i-1, ], nsubsamples)
    v <- c*v + (1-c)*g^2
    theta[i, ] <- theta[i-1, ] - a*g / (sqrt(v + epsilon))
    cost[i] <- f(y.samples, X.samples, theta[i, ], nsubsamples)
  }
  return(list(theta=theta, cost=cost))
}


### Algorithm 7: Adam given a cost function f and its gradient 
adam <- function(f, grad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, c, epsilon, m0, v0) {
  theta <- matrix(data=NA, nrow=niters, ncol=npars)
  cost <- vector(mode="numeric", length=niters)
  m <- vector(mode="numeric", length=niters)
  v <- vector(mode="numeric", length=niters)

  # Starting values
  theta[1, ]<- theta0
  cost[1] <- f(y, X, theta0, ndata)
  m <- m0
  v <- v0 
  
  # Updating steps 
  for(i in 2:niters){
    dataset <- cbind(y, X)
    random.index <- dataset[sample(nrow(dataset), nsubsamples), ]
    X.samples <- cbind(rep(1, times=nsubsamples), random.index[ ,2:npars+1])
    y.samples <- random.index[ ,1]
    g <- grad(y.samples, X.samples, theta[i-1, ], nsubsamples)
    m <- b*m + (1-b)*g
    v <- c*v + (1-c)*g^2
    m.hat <- m / (1-(b^i))
    v.hat <- v / (1-(c^i))
    theta[i, ] <- theta[i-1, ] - a*m.hat / (sqrt(v.hat) + epsilon)
    cost[i] <- f(y.samples, X.samples, theta[i, ], nsubsamples)
  }
  return(list(theta=theta, cost=cost))
}
