### Load R scripts with optimizers and cost functions

source("optimizers.r")
source("cost.r")
library(ggplot2)

### Load data, standardize covariate and create design matrix

### Load CSV file
weather <- read.csv("weather_data.csv", header=TRUE)

### Subset data, only keep variables of interest   
weather <- cbind(weather[, c("apparent_temperature", "temperature", "humidity", 
                             "wind_speed", "wind_bearing", "pressure")])

### Standardisation of covariates 
weather[, 2:6] <- scale(weather[, 2:6])

### Number of data points 
ndata <- nrow(weather)

## Number of linear regression parameters, including intercept 
npars <- 6

## Create design matrix X, including a first column of ones for intercept
X <- as.matrix(cbind(rep(1, ndata), weather[, 2:6]))
y <- weather$apparent_temperature


### Fit multiple linear regression model to the data
lm_out <- lm(y ~., data=as.data.frame(X[, 2:6]))
summary(lm_out)

## Estimated Linear Regressiion Coefficients 
lm_out$coefficients


################################ RUN OPTIMIZING ALGORITHMS ############################################### 


### Initial values of parameters, common across optimization algorithms
theta0 <- c(-5, -3, 4, 1, 10, -9)


### ALGORITHM 1 - GD: GRADIENT DESCENT ###

# Number of data points
ndata <- nrow(weather)

# learning rate
a <- 0.6

# number of iterations
niters <- 70

# Run GD
gd_out <- gd(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters)



### ALGORITHM 2 - SGD: STOCHASTIC GRADIENT DESCENT ###

# learning rate
a <- 0.5

# number of iterations
niters <- 200

# number of subsamples per iteration
nsubsamples <- 1000

# Run SGD
sgd_out <- sgd(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters, nsubsamples)



### ALGORITHM 3 - MSGD: STOCHASTIC GRADIENT DESCENT MOMENTUM ###

# learning rate
a <- 0.26

# number of iterations
niters <- 200

# number of subsamples per iteration
nsubsamples <- 1000

# memory factor for "momentum"
b <- 0.75

# initial value of "momentum"
m0 <- rep(0, npars)

# Run MSGD
msgd_out <- msgd(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, m0)



### ALGORITHM 4 - NAGSGD: STOCHASTIC GRADIENT DESCENT NESTEROV ACCELERATED GRADIENT 

# learning rate
a <-  0.5

# number of iterations
niters <- 200

# number of subsamples per iteration
nsubsamples <- 1000

# memory factor for "momentum"
b <- 0.6

# initial value of "momentum"
m0 <- rep(0, npars)

# Run NAGSGD
nagsgd_out <- nagsgd(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters, nsubsamples, b, m0)



### ALGORITHM 5 - AdaGad: ADAPTIVE GRADIENT ALGORITHM ###

# learning rate
a <- 2.6

# number of iterations
niters <- 200

# number of subsamples per iteration
nsubsamples <- 1000

# smoothing term
epsilon <- 1e-8

# initial value of G
G0 <- rep(0, npars)

# Run AdaGrad
adagrad_out <- adagrad(lmf, lmgrad, y, X, theta0, npars, ndata, a, niters, nsubsamples, epsilon, G0)



### ALGORITHM 6 - RMSProp: ROOT MEAN SQUARE PROPAGATION ###

# learning rate
a <- 0.208

# number of iterations
niters <- 200

# number of subsamples per iteration
nsubsamples <- 1000

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
niters <- 200

# number of subsamples per iteration
nsubsamples <- 1000

# memory factor for "momentum"
b <- 0.6

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



#################################### Answwers to Questions 2a) - 2e) #######################################################


### QUESTION 2a): save design matrix X to file "answer2a.csv"

write.table(X, file="answer2a.csv", row.names=FALSE, col.names=FALSE, sep=",")


### QUESTION 2b): save parameter estimates to file "answer2b.csv"

theta_estimates <- cbind(
  LM = as.vector(lm_out$coefficients),
  GD = as.vector(tail(gd_out$theta, n=1)),
  SGD = as.vector(tail(sgd_out$theta, n=1)),
  MSGD = as.vector(tail(msgd_out$theta, n=1)),
  NAGSGD = as.vector(tail(nagsgd_out$theta, n=1)),
  AdaGrad = as.vector(tail(adagrad_out$theta, n=1)),
  RMSProp = as.vector(tail(rmsprop_out$theta, n=1)),
  Adam = as.vector(tail(adam_out$theta, n=1))
)
write.table(theta_estimates, file="answer2b.csv", row.names=FALSE, col.names=TRUE, sep=",", quote=FALSE)


### QUESTION 2c): save last values of cost function to file "answer2c.csv"

final_cost_values <- cbind(
  GD = tail(gd_out$cost, n=1),
  SGD = tail(sgd_out$cost, n=1),
  MSGD = tail(msgd_out$cost, n=1),
  NAGSGD = tail(nagsgd_out$cost, n=1),
  AdaGrad = tail(adagrad_out$cost, n=1),
  RMSProp = tail(rmsprop_out$cost, n=1),
  Adam = as.vector(tail(adam_out$cost, n=1))
)
write.table(final_cost_values, file="answer2c.csv", row.names=FALSE, col.names=TRUE, sep=",", quote=FALSE)


### QUESTION 2d): save plot of cost function to file "answer2d.pdf"

pdf(file="answer2d.pdf", height=8, width=10)

ggplot() +
  geom_path(data = as.data.frame(nagsgd_out$cost), aes(x = 1:niters, y = nagsgd_out$cost, color = "red")) +
  geom_path(data = as.data.frame(adagrad_out$cost), aes(x = 1:niters, y = adagrad_out$cost, color = "blue")) +
  ylim(0, 350) +
  labs(title="Q.2.d: Plot of cost function", x="Iterations", y="Cost", 
       subtitle="NAGSGD: a=0.5 , b=0.6  |  AdaGrad: a=2.6") +
  scale_color_identity(name = "Optimizer", breaks = c("red", "blue"), labels = c("NAGSGD", "AdaGrad"), guide = "legend") +
  theme(legend.title = element_text(size=9.4, face="bold")) 

dev.off()


### QUESTION 2e): save phase plot of parameter theta_2 vs parameter theta_3 to file "answer2e.pdf"

pdf(file="answer2e.pdf", height=8, width=10)

ggplot() +
  geom_path(data = as.data.frame(gd_out$theta), aes(x = gd_out$theta[, 4], y = gd_out$theta[, 3], color = "red")) +
  geom_path(data = as.data.frame(msgd_out$theta), aes(x = msgd_out$theta[, 4], y = msgd_out$theta[, 3], color = "blue")) +
  geom_path(data = as.data.frame(adagrad_out$theta), aes(x = adagrad_out$theta[, 4], y = adagrad_out$theta[, 3], color = "black")) +
  geom_point(aes(x=lm_out$coefficients[4], y=lm_out$coefficients[3]), size = 1.5, colour="green") +
  labs(title="Q.2.e: Phase plot of" ~ theta[2]~vs~theta[3], x=expression(theta[3]), y=expression(theta[2]), 
       subtitle="GD: a=0.6  |  MSGD: a=0.26 , b=0.75  |  AdaGrad: a=2.6") +
  scale_color_identity(name = "Optimizer", breaks = c("red", "blue", "black"), labels = c("GD", "MSGD", "AdaGrad"), guide = "legend") +
  theme(legend.title = element_text(size=9.4, face="bold")) 

dev.off()
