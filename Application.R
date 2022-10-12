# Application Study
library(mvtnorm)

source("/Users/kdavis07/Documents/GitHub/CoPrimarySWCRT/powerSampleCal_IU_HoopGir.R")
source("/Users/kdavis07/Documents/GitHub/CoPrimarySWCRT/powerSampleCal_omnibus_HoopGir.R")

# Power under IU-test
## Parameter Inputs
sd1<-sqrt(611.13)
sd2<-sqrt(695.73)
rho2<-matrix(c(1,0.58,0.58,1),2)
rho0<-matrix(c(0.006,0,0,0.029),2);rho1<-matrix(c(0.00002,0,0,0.0068),2);deltas<-c(0.3*sd1,0.35*sd2);t<-5;N=(t-1)*4;m<-12

calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
# I=16 and N=12 gives 86.3% power


##### Sensitivity Analysis for Application Study - IU-test #####
### CAC=0.2
rho1<-matrix(c(0.00122,0,0,0.0059),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### CAC=0.5
rho1<-matrix(c(0.003,0,0,0.015),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{12} = -0.002
rho1<-matrix(c(0.00122,0,0,0.0059),2) #reset
rho0<-matrix(c(0.006,-0.002,-0.002,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{12} = 0.002
rho0<-matrix(c(0.006,0.002,0.002,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{2} = 0.023
rho0<-matrix(c(0.006,0,0,0.023),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{2} = 0.035
rho0<-matrix(c(0.006,0,0,0.035),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.005
rho0<-matrix(c(0.005,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.007
rho0<-matrix(c(0.007,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.46
rho0<-matrix(c(0.006,0,0,0.029),2) #reset
rho2<-matrix(c(1,0.46,0.46,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.7
rho2<-matrix(c(1,0.7,0.7,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)

#Additional scenarios (up to 40%) - don't forget to reset parameters before running
### CAC=0
rho1<-matrix(c(0,0,0,0),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### CAC=0.8
rho1<-matrix(c(0.0049,0,0,0.0235),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{12} = -0.004
rho1<-matrix(c(0.00122,0,0,0.0059),2) #reset
rho0<-matrix(c(0.006,-0.004,-0.004,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{12} = 0.004
rho0<-matrix(c(0.006,0.004,0.004,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{2} = 0.017
rho0<-matrix(c(0.006,0,0,0.017),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{2} = 0.041
rho0<-matrix(c(0.006,0,0,0.041),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.004
rho0<-matrix(c(0.004,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.008
rho0<-matrix(c(0.008,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.35
rho0<-matrix(c(0.006,0,0,0.029),2) #reset
rho2<-matrix(c(1,0.35,0.35,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.81
rho2<-matrix(c(1,0.81,0.81,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)

#Additional scenarios x2 (up to 60%) - don't forget to reset parameters before running
### rho_0^{2} = 0.012
rho0<-matrix(c(0.006,0,0,0.012),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{2} = 0.046
rho0<-matrix(c(0.006,0,0,0.046),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.002
rho0<-matrix(c(0.002,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.010
rho0<-matrix(c(0.010,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.23
rho0<-matrix(c(0.006,0,0,0.029),2) #reset
rho2<-matrix(c(1,0.23,0.23,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.93
rho2<-matrix(c(1,0.93,0.93,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)



# Power under omnibus test for SW-CRT
## Reset Parameter Inputs
sd1<-sqrt(611.13)
sd2<-sqrt(695.73)
rho2<-matrix(c(1,0.58,0.58,1),2)
rho0<-matrix(c(0.006,0,0,0.029),2);rho1<-matrix(c(0.00002,0,0,0.0068),2);deltas<-c(0.052*sd1,0.102*sd2);t<-5;N=(t-1)*4;m<-12

calPower_omnibus(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho0,rho1,rho2,N,t,m,K=2,alpha=0.05)
# Can test for standardized effect sizes as small as 0.052 and 0.102 with 86.5% power.



#Comparison to parallel-arm CRT under IU-test (using Yang et al)
# Power function under IU-test for a parallel-arm CRT - taken by Yang et al.
calPower_ttestIU <- function(betas,deltas,vars,rho01,rho2,N,r,m,K,alpha)
{
  #variance of trt assignment
  sigmaz.square <- r*(1-r) 
  
  ####Define function to calculate correlation between betas#####
  calCovbetas <- function(vars,rho01,rho2, sigmaz.square, m, K){
    rho0k <- diag(rho01)
    sigmak.square <-(1+(m-1)*rho0k)*vars/(m*sigmaz.square)
    covMatrix <- diag(sigmak.square)
    for(row in 1:K )
    {
      for(col in 1:K)
      {
        if(row != col){
          covMatrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]+(m-1)*rho01[row,col])/(m*sigmaz.square)
        }
      }
    }
    return(covMatrix)  
  }
  
  
  ####Define function to calculate correlation between test statistics #####
  calCorWks <-  function(vars,rho01,rho2, sigmaz.square, m, K)
  {
    top <- calCovbetas(vars,rho01,rho2, sigmaz.square, m, K)
    wCor <- diag(K)
    for(row in 1:K )
    {
      for(col in 1:K)
      {
        if(row != col){
          wCor[row,col] <- top[row,col]/sqrt(top[row,row]*top[col,col])
        }
      }
    }
    return(wCor)
  }
  sigmaks.sq <- diag(calCovbetas(vars,rho01,rho2, sigmaz.square, m, K))
  meanVector <- sqrt(N)*(betas-deltas)/sqrt(sigmaks.sq)
  wCor <- calCorWks(vars,rho01,rho2, sigmaz.square, m, K)
  criticalValue <- qt(p=(1-alpha), df=(N-2*K))
  pred.power <- pmvt(lower = rep(criticalValue,K),upper=rep(Inf,K),df = (N-2*K) , sigma = wCor,delta=meanVector)[1]
  return(pred.power)
}

# What if we assume a parallel-arm CRT design in our study using our current inputs?
## Parameter Inputs
sd1<-sqrt(611.13)
sd2<-sqrt(695.73)
rho2<-matrix(c(1,0.58,0.58,1),2)
rho01<-matrix(c(0.006,0,0,0.029),2);betas<-c(0.3*sd1,0.35*sd2);t<-5;N=(t-1)*4;m<-12*t

calPower_ttestIU(betas,deltas=c(0,0),vars=c(sd1^2,sd2^2),rho01=rho01,rho2=rho2,N,r=0.5,m,K=2,alpha=0.05)
# Same design parameters in a parallel-arm CRT has 91.5% power under IU-test



#Comparison to parallel-arm CRT under omnibus test (using Yang et al)
# Power function under omnibus test for a parallel-arm CRT - taken by Yang et al.
constrRiE <- function(rho01,rho2,K,vars)
{ rho0k <- diag(rho01)
SigmaE_Matrix <- diag((1-rho0k)*vars)
for(row in 1:K )
{
  for(col in 1:K)
  {
    if(row != col){
      SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]-rho01[row,col])
    }
  }
}
return(SigmaE_Matrix)
}

#####function to construct covariance matrix Sigma_phi for Y_i########
constrRiP <- function(rho01,K,vars)
{ rho0k <- diag(rho01)
SigmaP_Matrix <- diag(rho0k*vars)
for(row in 1:K )
{
  for(col in 1:K)
  {
    if(row != col){
      SigmaP_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho01[row,col]
    }
  }
}
return(SigmaP_Matrix)
}
####Define function to calculate covariance between betas#####
calCovbetas <- function(vars,rho01,rho2,cv, sigmaz.square, m, K){
  sigmaE <- constrRiE(rho01,rho2,K,vars)
  sigmaP <- constrRiP(rho01,K,vars)
  tmp <- solve(diag(1,K)-cv^2*(m*sigmaP %*% solve(sigmaE + m*sigmaP) %*% sigmaE %*% solve(sigmaE + m*sigmaP) ))
  covMatrix <- 1/(m*sigmaz.square)*(sigmaE+m*sigmaP)%*%tmp
  covMatrix <- (covMatrix +t(covMatrix))/2  # symmerize the off-diagonal
  return(covMatrix)  
}

calPower_Omnibus <- function(beta, vars,rho01,rho2,N, cv, m,r=0.5, K=2, alpha=0.05){
  sigmaz.square <- r*(1-r) 
  omega <- calCovbetas(vars,rho01,rho2,cv, sigmaz.square, m, K)
  beta<-matrix(beta,2,1)
  #Variance Inflation of coprimary outcomes
  tau <- N*t(beta) %*% solve(omega) %*% beta
  Fscore <- qf(1-0.05, df1=K, df2=N-2*K, ncp=0, lower.tail = TRUE, log.p = FALSE)
  power <-1- pf(Fscore, df1=K, df2=N-2*K, ncp=tau, lower.tail = TRUE, log.p = FALSE)
  return(power)
}

# What if we assume a parallel-arm CRT design in our study using our current inputs for omnibus?
## Parameter Inputs
sd1<-sqrt(611.13)
sd2<-sqrt(695.73)
rho2<-matrix(c(1,0.58,0.58,1),2)
rho01<-matrix(c(0.006,0,0,0.029),2);beta<-c(0.052*sd1,0.102*sd2);t<-5;N=(t-1)*4;m<-12*t

calPower_Omnibus(beta,vars=c(sd1^2,sd2^2),rho01=rho01,rho2=rho2,N,cv=0,m,r=0.5,K=2,alpha=0.05)
# Using standardized effect sizes of 0.052 and 0.102 gives only 12.0% power compared to the 86.5% power under SW-CRT.





