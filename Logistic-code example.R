# original link below ####
### https://tomroth.com.au/logistic/

# simulate the data ####
set.seed(2016)
#simulate data 
#independent variables
x1 = rnorm(30,3,2) + 0.1*c(1:30)
x2 = rbinom(30, 1,0.3)
x3 = rpois(n = 30, lambda = 4)
x3[16:30] = x3[16:30] - rpois(n = 15, lambda = 2)

#dependent variable 
y = c(rbinom(5, 1,0.1),rbinom(10, 1,0.25),rbinom(10, 1,0.75),rbinom(5, 1,0.9))

# design matrix ####
x0 = rep(1,30) #bias
X = cbind(x0,x1,x2,x3)

# logistic code ####
manual_logistic_regression = function(X,y,threshold = 1e-10, max_iter = 100)
  #A function to find logistic regression coeffiecients 
  #Takes three inputs: 
{
  #A function to return p, given X and beta
  #We'll need this function in the iterative section
  calc_p = function(X,beta)
  {
    beta = as.vector(beta)
    return(exp(X%*%beta) / (1+ exp(X%*%beta)))
  }  
  
  #### setup bit ####
  
  #initial guess for beta
  beta = rep(0,ncol(X))
  
  #initial value bigger than threshold so that we can enter our while loop 
  diff = 10000 
  
  #counter to ensure we're not stuck in an infinite loop
  iter_count = 0
  
  #### iterative bit ####
  while(diff > threshold ) #tests for convergence
  {
    #calculate probabilities using current estimate of beta
    p = as.vector(calc_p(X,beta))
    
    #calculate matrix of weights W
    W =  diag(p*(1-p)) 
    
    #calculate the change in beta
    beta_change = solve(t(X)%*%W%*%X) %*% t(X)%*%(y - p)
    
    #update beta
    beta = beta + beta_change
    
    #calculate how much we changed beta by in this iteration 
    #if this is less than threshold, we'll break the while loop 
    diff = sum(beta_change^2)
    
    #see if we've hit the maximum number of iterations
    iter_count = iter_count + 1
    if(iter_count > max_iter) {
      stop("This isn't converging, mate.")
    }
  }
  #make it pretty 
  coef = c("(Intercept)" = beta[1], x1 = beta[2], x2 = beta[3], x3 = beta[4])
  return(coef)
}

#using R package 
M1 = glm(y~x1+x2+x3, family = "binomial")
M1$coefficients
# (Intercept)          x1          x2          x3 
# -1.3512086   0.3191309   0.2033449  -0.0832102 

#our version
manual_logistic_regression(X,y)
# (Intercept)          x1          x2          x3 
# -1.3512086   0.3191309   0.2033449  -0.0832102 

#### CODE2- via optim ####
# ===
# https://stats.stackexchange.com/questions/81000/calculate-coefficients-in-a-logistic-regression-with-r
# ===

# compute the logistic regression parameters as 
#   an optimal value
#==
# define the logistic transformation
logit = function(mX, vBeta) {
  return(exp(mX %*% vBeta)/(1+ exp(mX %*% vBeta)) )
}

# stable parametrisation of the log-likelihood function
# Note: The negative of the log-likelihood is being returned, since we will be
# /minimising/ the function.
logLikelihoodLogitStable = function(vBeta, mX, vY) {
  return(-sum(
    vY*(mX %*% vBeta - log(1+exp(mX %*% vBeta))) + (1-vY)*(-log(1 + exp(mX %*% vBeta)))
  )) 
}

# initial set of parameters
# tbeta = c(10, -0.1, -0.3, 0.001) # arbitrary starting parameters

# minimise the (negative) log-likelihood to get the logit fit
optimLogit = optim(tbeta, logLikelihoodLogitStable, mX = X, vY = y, method = 'BFGS', hessian=TRUE)
logit(X, optimLogit$par)

# LEo-CODE ####
I <- nrow(X)
y.I <- array(y, I)
X.I <- aperm(X, c(2,1))
f <- function(theta, X.I, y.I){
  # prior<- sum(log(dnorm(theta))) #using the prior throws the estimated elsewhere
  theta.I <- theta%o%rep(1,I)
  ll.I<-y.I*(apply(X.I*theta.I,2,sum) - log(1+ exp(apply(X.I*theta.I,2,sum)))) + 
    (1-y.I)*(-log(1+ exp(apply(X.I*theta.I,2,sum))))
  # f<--2*(prior+sum(ll.I))
  f<- -sum(ll.I)
}

m.f <- optim(tbeta, f, X.I = X.I, y.I = y.I, method = 'BFGS', hessian=TRUE)
m.f$par
optimLogit$par
M1$coefficients
manual_logistic_regression(X,y)

#### another data ####
mydata <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
mylogit <- glm(admit ~ gre + gpa, data = mydata, family = "binomial")
summary(mylogit)

X2 <- as.matrix(cbind(1, mydata[,2:3]))
y2 <- as.matrix(mydata[,1])


manual_logistic_regression(X2, y2)


### REFINED CODE1 ####

#A function to return p, given X and beta
#We'll need this function in the iterative section
calc_p = function(X,beta)
{
  beta = as.vector(beta)
  return(exp(X%*%beta) / (1+ exp(X%*%beta)))
}  

#initial guess for beta
beta = function(X){rep(0,ncol(X))}

log_reg = function(X,y, beta, threshold = 1e-10, max_iter = 100)
  #Takes three inputs: 
{
  #### setup bit ####
  #initial value bigger than threshold so that we can enter our while loop 
  diff = 10000 
  
  #counter to ensure we're not stuck in an infinite loop
  iter_count = 0
  
  #### iterative bit ####
  while(diff > threshold ) #tests for convergence
  {
    #calculate probabilities using current estimate of beta
    p = as.vector(calc_p(X,beta))
    
    #calculate matrix of weights W
    W =  diag(p*(1-p)) 
    
    #calculate the change in beta
    beta_change = solve(t(X)%*%W%*%X) %*% t(X)%*%(y - p)
    
    #update beta
    beta = beta + beta_change
    
    #calculate how much we changed beta by in this iteration 
    #if this is less than threshold, we'll break the while loop 
    diff = sum(beta_change^2)
    
    #see if we've hit the maximum number of iterations
    iter_count = iter_count + 1
    if(iter_count > max_iter) {
      stop("This isn't converging, mate.")
    }
  }
  return(beta)
}


log_reg(X=X2, y=y2, beta= beta(X2))
coef(mylogit)
tbeta2 <-c(-20, 0.5,0.6)
optimLogit2 = optim(tbeta2, logLikelihoodLogitStable, mX = X2, vY = y2, method = 'BFGS', hessian=TRUE)
optimLogit2$par
logit(X2, optimLogit2$par)
