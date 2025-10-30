##############################################
##############################################
#### Tri-variate Chi-square CDF ##############
##############################################
##############################################

library(MASS)
library(matrixcalc)

#### generate normal vector
N = 100000
mu=rep(0,5*3)

propmat <- matrix(c(1, sqrt(1/2), sqrt(1/3), sqrt(1/2), 1, sqrt(2/3), sqrt(1/3), sqrt(2/3), 1), nrow=3)
I5 <- diag(rep(1,5), nrow=5)
Sigma <- direct.prod(propmat, I5)

z <- mvrnorm(n=N, mu=mu, Sigma=Sigma)
z1 <- z[,1:5]
z2 <- z[,6:10]
z3 <- z[,10:15]

#### transform to quadratic form
x1 <- rep(NA, N)
x2 <- rep(NA, N)
x3 <- rep(NA, N)
for (i in 1:N){
  x1[i] = t(z1[i,])%*%z1[i,]
  x2[i] = t(z2[i,])%*%z2[i,]
  x3[i] = t(z3[i,])%*%z3[i,]
}

empirical_cdf <- function(x1, x2, x3, c1, c2, c3, N){
  sub = rep(NA, N)
  
  for (i in 1:N){
    sub[i]=(x1[i]<=c1)*(x2[i]<=c2)*(x3[i]<=c3)
  }
  
  a = sum(sub)/N
  return(a)
}

x_search <- function(start, end, by, x1, x2, x3, alpha, OBF=F){
  x=seq(start,end,by)
  y = sapply(x, function(x){empirical_cdf(x1=x1, x2=x2, x3=x3, c1=x*ifelse(OBF==F,1,sqrt(3)), c2=x*ifelse(OBF==F,1,sqrt(3/2)), c3=x, N=N)})
  index=which(y>(1-alpha))[1]
  
  if(OBF==F){c=rep(x[index],3)}else{c=c(x[index]*sqrt(3),x[index]*sqrt(2),x[index])}
  
  return(c)
}

##
x_search(start=10, end=30, by=0.01, x1=x1, x2=x2, x3=x3, alpha=0.05, OBF=F) 
x_search(start=10, end=30, by=0.01, x1=x1, x2=x2, x3=x3, alpha=0.05, OBF=T) 
