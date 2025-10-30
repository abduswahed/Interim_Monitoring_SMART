##############################################
##############################################
#### Example code to analyze temp.RData ######
##############################################
library('MASS')
load('temp.RData')

###### observed treatment means
mua1b1.hat = mean(temp[temp$trt=='A1B1',]$y)
mua1b2.hat = mean(temp[temp$trt=='A1B2',]$y)
mua2b1.hat = mean(temp[temp$trt=='A2B1',]$y)
mua2b2.hat = mean(temp[temp$trt=='A2B2',]$y)
mua1c1.hat = mean(temp[temp$trt=='A1C1',]$y)
mua1c2.hat = mean(temp[temp$trt=='A1C2',]$y)
mua2c1.hat = mean(temp[temp$trt=='A2C1',]$y)
mua2c2.hat = mean(temp[temp$trt=='A2C2',]$y)
###### observed treatment variance
sda1b1.hat = sd(temp[temp$trt=='A1B1',]$y)
sda1b2.hat = sd(temp[temp$trt=='A1B2',]$y)
sda2b1.hat = sd(temp[temp$trt=='A2B1',]$y)
sda2b2.hat = sd(temp[temp$trt=='A2B2',]$y)
sda1c1.hat = sd(temp[temp$trt=='A1C1',]$y)
sda1c2.hat = sd(temp[temp$trt=='A1C2',]$y)
sda2c1.hat = sd(temp[temp$trt=='A2C1',]$y)
sda2c2.hat = sd(temp[temp$trt=='A2C2',]$y)

###### observed randomization probabilities
prob1.hat = 1-sum(temp$A2==1)/dim(temp)[1] ## probability to A1
p1.hat = 1-sum(temp$B2==1, na.rm=T)/sum(temp$R==1) ## probability to B1 for responder
q1.hat = 1-sum(temp$C2==1, na.rm=T)/sum(temp$R==0) ## probability to C1 for non-responder

kappa=c(1/prob1.hat, 1/(1-prob1.hat))
p=c(p1.hat, 1-p1.hat)
q=c(q1.hat, 1-q1.hat)

##### first stage response rate
pi1.hat = sum(temp[temp$A2==0,]$R==1)/sum(temp$A2==0)
pi2.hat = sum(temp[temp$A2==1,]$R==1)/sum(temp$A2==1)
pi.hat=c(pi1.hat,pi2.hat)

##### strategy means
mean_array = array(c(mua1b1.hat, mua1b2.hat, mua2b1.hat, mua2b2.hat, mua1c1.hat, mua1c2.hat, mua2c1.hat, mua2c2.hat), dim=c(2,2,2))

for (i in 1:2){
  for (j in 1:2){
    for (k in 1:2){
      ## Equation (3)
      W <- (temp$A2==(i-1))*{temp$R*(temp$B2==(j-1))/p[j]+(1-temp$R)*(temp$C2==(k-1))/q[k]}
      assign(paste("mu",i,j,k, sep=""), sum(W*temp$y)/sum(W))
    }
  }
}
mu_array = array(c(mu111, mu211, mu121, mu221, mu112, mu212, mu122, mu222), dim=c(2,2,2))

##### strategy variance
sd_array = array(c(sda1b1.hat, sda1b2.hat, sda2b1.hat, sda2b2.hat, sda1c1.hat, sda1c2.hat, sda2c1.hat, sda2c2.hat), dim=c(2,2,2))

for (i in 1:2){
  for (j in 1:2){
    for (k in 1:2){
      ## Equation (6)              
      W <-{temp$R*(temp$B2==(j-1))/p[j]+(1-temp$R)*(temp$C2==(k-1))/q[k]}
      W <- W[(temp$A2==(i-1))]
      y <- temp$y[(temp$A2==(i-1))]
      assign(paste("sigma", i, j, k, sep=""),  sum(W^2*(y-mu_array[i,j,k])^2)/sum(temp$A2==(i-1))/(sum(temp$A2==(i-1))-1))
    }
  }
}

for (i in 1:2){
  for (j in 1:2){
    ## Equation (7)
    W1 <- temp$R*(temp$B2==(j-1))
    W1 <- W1[(temp$A2==(i-1))]
    y <- temp$y[(temp$A2==(i-1))]
    assign(paste("sigma", i,j,1,"_",i,j,2, sep=""),
           sum(W1*(y-mu_array[i,j,1])*(y-mu_array[i,j,2]))/p[j]^2/sum(temp$A2==(i-1))^2)
    
    W3 <- (1-temp$R)*(temp$C2==(j-1))
    W3 <- W3[(temp$A2==(i-1))]
    assign(paste("sigma",i,1,j,"_",i,2,j, sep=""),
           sum(W3*(y-mu_array[i,1,j])*(y-mu_array[i,2,j]))/q[j]^2/sum(temp$A2==(i-1))^2)
  }
}

Sigma = matrix(c(sigma111, sigma111_112, sigma111_121, rep(0,5), 
                 sigma111_112, sigma112, 0, sigma112_122, rep(0,4),
                 sigma111_121, 0, sigma121, sigma121_122, rep(0,5), 
                 sigma112_122, sigma121_122, sigma122, rep(0,4),
                 rep(0,4), sigma211, sigma211_212, sigma211_221, 
                 rep(0,5), sigma211_212, sigma212, 0, sigma212_222, 
                 rep(0,4), sigma211_221, 0, sigma221, sigma221_222, 
                 rep(0,5), sigma212_222, sigma221_222, sigma222), nrow=8)

# all strategy means
mu = c(mu111, mu112, mu121, mu122, mu211, mu212, mu221, mu222)

# variance matrix Sigma adjusted for the inflation factor
k=21 # total number of parameters estimated
IF=dim(temp)[1]/(dim(temp)[1]-k)
Sigma = Sigma*IF

# Contrast matrix
C = matrix(c(rep(1,7), rep(c(-1, rep(0,7)), 6), -1), nrow=7)

# Wald-type test statistic in Equation (8)
t(C%*%mu)%*%ginv(C%*%Sigma%*%t(C))%*%C%*%mu
