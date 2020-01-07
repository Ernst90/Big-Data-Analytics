## Big Data Analytics (STATS 5016) - Masters Level 
## Project: Stochastic Optimization
## GUID: 2383746W

### Load R scripts with optimizers and cost functions

source("optimizers.r")
source("cost.r")
library(ggplot2)

### Load data, standardize covariate and create design matrix

## Load CSV file
simulated_data <- read.csv("simulated_data.csv", header=TRUE)

## Standardize covariate
simulated_data$covariate <- scale(simulated_data$covariate)

## Number of data points
ndata <- nrow(simulated_data)

## Number of linear regression parameters, including intercept 
npars <- 2

## Create design matrix X, including a first column of ones for intercept
X <- cbind(rep(1, times=ndata), simulated_data$covariate)
y <- simulated_data$y

### Fit simple linear regression model to the data
lm_out <- lm(y~., data=simulated_data)
summary(lm_out)

## Estimated Linear Regressiion Coefficients 
lm_out$coefficients


################################ RUN OPTIMIZING ALGORITHMS ############################################### 

## Initial values of parameters, common across optimization algorithms
theta0 <- c(7, -8)


### ALGORITHM 1 - GD: GRADIENT DESCENT ###

# learning rate
a <- 0.28

# number of iterations
niters <- 100

# Run GD
gd_out <- gd(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters)


### ALGORITHM 2 - SGD: STOCHASTIC GRADIENT DESCENT ###

# learning rate
a <- 0.085

# number of iterations
niters <- 100

# number of subsamples per iteration
nsubsamples <- 100

# Run SGD
sgd_out <- sgd(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters, nsubsamples)


### ALGORITHM 3 - MSGD: STOCHASTIC GRADIENT DESCENT MOMENTUM ###

# learning rate
a <- 0.105

# number of iterations
niters <- 100

# number of subsamples per iteration
nsubsamples <- 100

# memory factor for "momentum"
b <- 0.27

# initial value of "momentum"
m0 <- rep(0, npars)

# Run MSGD
msgd_out <- msgd(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, m0)


### ALGORITHM 4 - NAGSGD: STOCHASTIC GRADIENT DESCENT NESTEROV ACCELERATED GRADIENT 

# learning rate
a <- 0.155

# number of iterations
niters <- 100

# number of subsamples per iteration
nsubsamples <- 100

# memory factor for "momentum"
b <- 0.58

# initial value of "momentum"
m0 <- rep(0, npars)

# Run NAGSGD
nagsgd_out <- nagsgd(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, m0)


### ALGORITHM 5 - AdaGad: ADAPTIVE GRADIENT ALGORITHM ###

# learning rate
a <- 1.5

# number of iterations
niters <- 100

# number of subsamples per iteration
nsubsamples <- 100

# smoothing term
epsilon <- 1e-8

# initial value of G
G0 <- rep(0, npars)

# Run AdaGrad
adagrad_out <- adagrad(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters, nsubsamples, epsilon, G0)


### ALGORITHM 6 - RMSProp: ROOT MEAN SQUARE PROPAGATION ###

# learning rate
a <- 0.054

# number of iterations
niters <- 100

# number of subsamples per iteration
nsubsamples <- 100

# memory factor for squared gradient
c <- 0.999

# smoothing term
epsilon <- 1e-8

# initial value of squared gradient
v0 <- rep(0, npars)

# Run RMSprop
rmsprop_out <- rmsprop(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters, nsubsamples, c, epsilon, v0)


### ALGORITHM 7 - ADAM: ADAPTIVE MOMENT ESTIMATION ###

# learning rate
a <- 0.45

# number of iterations
niters <- 100

# number of subsamples per iteration
nsubsamples <- 100

# memory factor for "momentum"
b <- 0.554

# memory factor for squared gradient
c <- 0.999

# smoothing term
epsilon <- 1e-8

# initial value of "momentum"
m0 <- rep(0, npars)

# initial value of squared gradient
v0 <- rep(0, npars)

# Run Adam
adam_out <- adam(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, c, epsilon, m0, v0)


#################################### Answwers to Questions 1a) - 1e) #######################################################

### QUESTION 1a): save design matrix X to file "answer1a.csv"

write.table(X, file="answer1a.csv", row.names=FALSE, col.names=FALSE, sep=",")

### QUESTION 1b): save parameter estimates to file "answer1b.csv"

theta_estimates <- cbind(
 LM = as.vector(lm_out$coefficients),
 GD = as.vector(tail(gd_out$theta, n=1)),
 SGD = as.vector(tail(sgd_out$theta, n=1)),
 MSGD = as.vector(tail(msgd_out$theta, n=1)),
 NAGSGD = as.vector(tail(nagsgd_out$theta, n=1)),
 AdaGrad = as.vector(tail(adagrad_out$theta, n=1)),
 RMSProp = as.vector(tail(rmsprop_out$theta, n=1)),
 ADAM = as.vector(tail(adam_out$theta, n=1))
)

write.table(theta_estimates, file="answer1b.csv", row.names=FALSE, col.names=TRUE, sep=",", quote=FALSE)

### QUESTION 1c): save last values of cost function to file "answer1c.csv"

final_cost_values <- cbind(
 GD = tail(gd_out$cost, n=1),
 SGD = tail(sgd_out$cost, n=1),
 MSGD = tail(msgd_out$cost, n=1),
 NAGSGD = tail(nagsgd_out$cost, n=1),
 AdaGrad = tail(adagrad_out$cost, n=1),
 RMSProp = tail(rmsprop_out$cost, n=1),
 ADAM = tail(adam_out$cost, n=1)
)
write.table(final_cost_values, file="answer1c.csv", row.names=FALSE, col.names=TRUE, sep=",", quote=FALSE)

### QUESTION 1d): save plot of cost function to file "answer1d.pdf"

pdf(file="answer1d.pdf", height=8, width=10)

# NAGSGD vs AdaGrad
ggplot() +
  geom_path(data = as.data.frame(nagsgd_out$cost), aes(x = 1:niters, y = nagsgd_out$cost, color = "orange")) +
  geom_path(data = as.data.frame(adagrad_out$cost), aes(x = 1:niters, y = adagrad_out$cost, color = "purple")) +
  xlim(0, 75) + ylim(0, 80) + theme_bw() +
  labs(title="Q.1.d: Plot of cost function", x="Iterations", y="Cost", 
       subtitle="NAGSGD: a=0.155 , b=0.58  |  AdaGrad: a=1.5") +
  scale_color_identity(name = "Optimizer", breaks = c("orange", "purple"), 
                       labels = c("NAGSGD", "AdaGrad"), 
                       guide = "legend") +
  theme(legend.title = element_text(size=9.4, face="bold")) 

########### Additional 1d) - all other Algorithms ############
ggplot() +
  geom_path(data = as.data.frame(gd_out$cost), aes(x = 1:niters, y = gd_out$cost, color = "red")) +
  geom_path(data = as.data.frame(sgd_out$cost), aes(x = 1:niters, y = sgd_out$cost, color = "cyan")) +
  geom_path(data = as.data.frame(msgd_out$cost), aes(x = 1:niters, y = msgd_out$cost, color = "blue")) +
  geom_path(data = as.data.frame(rmsprop_out$cost), aes(x = 1:niters, y = rmsprop_out$cost, color = "brown")) +
  geom_path(data = as.data.frame(adam_out$cost), aes(x = 1:niters, y = adam_out$cost, color = "forestgreen")) +
  xlim(0, 75) + ylim(0, 80) + theme_bw() +
  labs(title="Q.1.d: Additional plots of cost functions", x="Iterations", y="Cost", 
       subtitle="GD: a=0.28  |  SGD: a=0.085  |  MSGD: a=0.105 , b=0.27     
RMSProp: a=0.054 , c=0.999  |  ADAM: a=0.45 , b=0.554 , c=0.999") +
  scale_color_identity(name = "Optimizer", breaks = c("red", "cyan", "blue", "brown", "forestgreen"), 
                       labels = c("GD", "SGD", "MSGD", "RMSProp", "ADAM"), 
                       guide = "legend") +
  theme(legend.title = element_text(size=9.4, face="bold")) 

dev.off()


### QUESTION 1e): save phase plot of parameter theta_0 vs parameter theta_1 to file "answer1e.pdf"

pdf(file="answer1e.pdf", height=8, width=10)

## GD vs RMSprop
ggplot() +
  geom_path(data = as.data.frame(gd_out$theta), aes(x = gd_out$theta[, 2], y = gd_out$theta[, 1], color = "red")) +
  geom_path(data = as.data.frame(rmsprop_out$theta), aes(x = rmsprop_out$theta[, 2], y = rmsprop_out$theta[, 1], color = "brown")) +
  geom_point(aes(x=lm_out$coefficients[2], y=lm_out$coefficients[1]), size = 2.5, colour="black") + theme_bw() +
  labs(title="Q.1.e: Phase plot of" ~ theta[0]~vs~theta[1], x=expression(theta[1]), y=expression(theta[0]), 
       subtitle="GD: a=0.28  |  RMSProp: a=0.054 , c=0.999") +
  scale_color_identity(name = "Optimizer", breaks = c("red", "brown"), 
                       labels = c("GD", "RMSProp"), 
                       guide = "legend") +
  theme(legend.title = element_text(size=9.4, face="bold")) 

########### Additional 1e) - all other Algorithms ############
ggplot() +
  geom_path(data = as.data.frame(sgd_out$theta), aes(x = sgd_out$theta[, 2], y = sgd_out$theta[, 1], color = "cyan")) +
  geom_path(data = as.data.frame(msgd_out$theta), aes(x = msgd_out$theta[, 2], y = msgd_out$theta[, 1], color = "blue")) +
  geom_path(data = as.data.frame(nagsgd_out$theta), aes(x = nagsgd_out$theta[, 2], y = nagsgd_out$theta[, 1], color = "orange")) +
  geom_path(data = as.data.frame(adagrad_out$theta), aes(x = adagrad_out$theta[, 2], y = adagrad_out$theta[, 1], color = "purple")) +
  geom_path(data = as.data.frame(adam_out$theta), aes(x = adam_out$theta[, 2], y = adam_out$theta[, 1], color = "forestgreen")) +
  geom_point(aes(x=lm_out$coefficients[2], y=lm_out$coefficients[1]), size = 2.5, colour="black") + theme_bw() +
  labs(title="Q.1.e: Additional phase plots of" ~ theta[0]~vs~theta[1], x=expression(theta[1]), y=expression(theta[0]), 
       subtitle="SGD: a=0.085  |  MSGD: a=0.105 , b=0.27  |  NAGSGD: a=0.155 , b=0.58    
AdaGrad: a=1.5  |  ADAM: a=0.45 , b=0.554 , c=0.999") +
  scale_color_identity(name = "Optimizer", breaks = c("cyan", "blue", "orange", "purple", "forestgreen"), 
                       labels = c("SGD", "MSGD", "NAGSGD", "AdaGrad", "ADAM"), 
                       guide = "legend") +
  theme(legend.title = element_text(size=9.4, face="bold")) 

dev.off()
