##############################################
##############################################
#### Bi-variate Chi-square CDF ###############
##############################################
##############################################
# Reference: Dickhaus, T. and Royen, T. (2015).  A survey on multivariate chi-square distributions andtheir applications in testing multiple hypotheses.Statistics49,427–454

Laguerre <- function(x, n, alpha){
  #### a function to calcualte Laguerre polynomial
  sum=0
  for (i in 0:n){
    sum=sum+(-x)^i*gamma(n+alpha+1)/gamma(n-i+1)/gamma(alpha+i+1)/gamma(i+1)
  }  
  return(sum)
}

Gamma_n <- function(x, n, alpha){
  #### a function to calculate nth derivative of Gamma cdf
  ifelse(n==0, pgamma(q=x, shape=alpha),
         gamma(alpha+1)*gamma(n)/gamma(alpha+n)*Laguerre(x=x, alpha=alpha, n=n-1)*dgamma(x=x, shape=alpha+1))
}

coef_a <- function(r, n){
  a <- rep(NA, length=n+1)
  a[1] = 1
  
  if (n > 0){
    for (i in 1:n){
      temp=0
      for (m in 0:(i-1)){
        temp = temp + a[i-m]*sum(r^(2*(m+1))) 
      }
      a[i+1] = 1/2/(i)*temp  
    }
  }
  
  return(a[n+1])
}

joint_cdf <- function(x1, x2, r, v, max.iter=100, delta=10^-6){
  c=0
  
  for (i in 0:max.iter){
    temp=c
    c = c + coef_a(r=r, n=i)*Gamma_n(x=x1/2, n=i, alpha=v/2)*Gamma_n(x=x2/2, n=i, alpha=v/2)
    
    if (abs(c-temp)<delta) break
  }
  
  result=list()
  result$c = c
  result$iteration = i
  return(result)
}

x_search <- function(start, end, by, R, alpha, K){
  x=seq(start,end,by)
  y = sapply(x, function(x){joint_cdf(x1=x*K, x2=x, r=sqrt(eigen(R%*%t(R))$values), v=dim(R)[1])$c})
  index=which(y>(1-alpha))[1]
  
  return(x[index])
}

##
proportion <- 0.5
R<-diag(rep(sqrt(proportion),5), nrow=5, ncol=5)
x_search(start=1, end=10, by=0.01, R=R, alpha=0.05, K=1)
x_search(start=1, end=10, by=0.01, R=R, alpha=0.05, K=sqrt(1/proportion))