library(mvtnorm)
########################################################################################################################################################
##Power/Sample Size Calculation based on the t test##
# INPUT
# deltas: (delta_1,...,delta_K), the vector of treatment effect for (1st,...,Kth) endpoints
# margins: (margin_1,...,margin_K), the vector of non-inferiority margins, when delta_1 = ... = delta_K = 0,
#          superiority tests are performed on all endpoints
# vars: (var_1,...,var_K), the vector of marginal variance for (1st,...,Kth) endpoints
# rho01: a K by K dimensional matrix for the correlation parameters (rho0^k) and (rho1^kk')
# For rho01: 
#           the diagonal elements correspond to rho0^k's 
#           the off-diagonal elements correspond to (rho1^kk')'s 
#           For example, rho01[1,1] corresponds to rho0^1, which is the ICC for the first endpoint
#                        rho01[1,2] corresponds to rho1^12, which is the correlation of outcomes between subjects on the 1st and 2nd endpoints    
# rho2: a K by K dimensional matrix for the correlation parameters (rho2^kk')
# For rho2:
#           the diagonal elements are 1
#           the off-diagonal elements correspond to (rho2^kk')'s
#           For example, rho2[1,2] corresponds to rho2^12, which is the correlation of outcomes within same subject on the 1st and 2nd endpoints    
# N: number of clusters
# t: number of time periods
# m: cluster size
# K: number of endpoints
# alpha: upper bound of type I error rates over the whole null space
########################################################################################################################################################
####Function to Calculate Power Given Design Configurations based on the t test and normal (intersection-union test)#######
####Critical values c_1,...,c_K are set to t_alpha, (1-alpha)th quantile of the t distribution with df = N-2K###
calPower_IU <- function(deltas,margins,vars,rho01,rho2,N,t,m,K,alpha)
{
  # Create X matrix
  X<-NULL
  trtSeq<-matrix(0,t-1,t)
  trtSeq[upper.tri(trtSeq)]<-1
  g<-N/(t-1) # number of clusters per step
  X=trtSeq[rep(1:nrow(trtSeq),each=g),]
  
  # constant calculation
  U=sum(X)
  V=sum((rowSums(X))^2)
  W=sum((colSums(X))^2)
  
  #####function to construct covariance matrix Sigma_E for Y_i########
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
  ####Define function to calculate covariance between deltas#####
  calCovbetas <- function(vars,rho01,rho2){
    sigmaE <- constrRiE(rho01,rho2,K,vars)
    sigmaP <- constrRiP(rho01,K,vars)
    #tmp <- solve(diag(1,K)-cv^2*(m*sigmaP %*% solve(sigmaE + m*sigmaP) %*% sigmaE %*% solve(sigmaE + m*sigmaP) ))
    #covMatrix <- 1/(m*sigmaz.square)*(sigmaE+m*sigmaP)%*%tmp
    #covMatrix <- (covMatrix +t(covMatrix))/2  # symmerize the off-diagonal
  #Stepped wedge
    covMatrix <- (N/(m*(N*U-W)))*(sigmaE-((U^2-N*V)/(U^2+N*t*U-t*W-N*V))*solve(solve(sigmaE)+((N*U-W)/(m*(U^2+N*t*U-t*W-N*V)))*solve(sigmaP)))
    #covMatrix <- N*t*solve(m*(U^2+N*t*U-t*W-N*V)*solve(sigmaE)+(N*V-U^2)*solve(t*sigmaP+(1/m)*sigmaE))
    return(covMatrix)  
  }
  
  
  ####Define function to calculate correlation between test statistics #####
  calCorWks <-  function(vars,rho01,rho2)
  {
    top <- calCovbetas(vars,rho01,rho2)
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
  sigmaks.sq <- diag( calCovbetas(vars,rho01,rho2))
  #meanVector <- sqrt(N)*(deltas-margins)/sqrt(sigmaks.sq)
  meanVector <- (deltas-margins)/sqrt(sigmaks.sq)
  wCor <- calCorWks(vars,rho01,rho2)
  criticalValue.t <- qt(p=(1-alpha), df=(N-2*K))
  criticalValue.n <- qnorm(p=(1-alpha))
  pred.power.t <- pmvt(lower = rep(criticalValue.t,K),upper=rep(Inf,K),df = (N-2*K), sigma = wCor,delta=meanVector)[1]
  pred.power.n <- pmvnorm(lower = rep(criticalValue.n,K),upper=rep(Inf,K), sigma = wCor,mean=meanVector)[1]
  
  param <- list(vard=c(sigmaks.sq),pred.power.t=pred.power.t,pred.power.z=pred.power.n)
  return(param)
}

#Prelim Study
#deltas<-c(0.8,0.3); margins<-c(0,0); vars<-c(1.1,1.05); rho01<-matrix(c(0.02,0.01,0.01,0.015),2); rho2<-matrix(c(1,0.05,0.05,1),2); t<-4; N=(t-1)*5; m<-15;K<-2; alpha<-0.05
#deltas<-c(0.25,1.3); margins<-c(0,0); vars<-c(1.5,1.2); rho01<-matrix(c(0.2,0.1,0.1,0.15),2); rho2<-matrix(c(1,0.3,0.3,1),2); t<-5; N=(t-1)*7; m<-12;K<-2; alpha<-0.05
#deltas<-c(1.1,0.95); margins<-c(0,0); vars<-c(0.8,1.5); rho01<-matrix(c(0.04,0.02,0.02,0.05),2); rho2<-matrix(c(1,0.09,0.09,1),2); t<-3; N=(t-1)*4; m<-10;K<-2; alpha<-0.05
#deltas<-c(2,0.8); margins<-c(0,0); vars<-c(2,1.5); rho01<-matrix(c(0.06,0.02,0.02,0.04),2); rho2<-matrix(c(1,0.1,0.1,1),2); t<-3; N=(t-1)*8; m<-5;K<-2; alpha<-0.05
#deltas<-c(0.75,0.75); margins<-c(0,0); vars<-c(1.1,0.7); rho01<-matrix(c(0.1,0.05,0.05,0.15),2); rho2<-matrix(c(1,0.2,0.2,1),2); t<-4; N=(t-1)*4; m<-4;K<-2; alpha<-0.05

#Sim Study Development
rho2<-matrix(c(1,0.2,0.2,1),2)
  rho01<-matrix(c(0.01,0.005,0.005,0.01),2);deltas<-c(0.42,0.28);t<-3;N=(t-1)*4;m<-82
  rho01<-matrix(c(0.01,0.005,0.005,0.05),2);deltas<-c(0.25,0.35);t<-5;N=(t-1)*2;m<-50
  rho01<-matrix(c(0.01,0.005,0.005,0.1),2);deltas<-c(0.22,0.50);t<-4;N=(t-1)*3;m<-68
  
  rho01<-matrix(c(0.05,0.005,0.005,0.01),2);deltas<-c(0.33,0.31);t<-5;N=(t-1)*2;m<-37
  rho01<-matrix(c(0.05,0.025,0.025,0.05),2);deltas<-c(0.47,0.34);t<-4;N=(t-1)*4;m<-20
  rho01<-matrix(c(0.05,0.025,0.025,0.1),2);deltas<-c(0.26,0.50);t<-3;N=(t-1)*6;m<-55
  
  rho01<-matrix(c(0.1,0.005,0.005,0.01),2);deltas<-c(0.45,0.1);t<-5;N=(t-1)*5;m<-85
  rho01<-matrix(c(0.1,0.025,0.025,0.05),2);deltas<-c(0.49,0.42);t<-3;N=(t-1)*5;m<-29
  rho01<-matrix(c(0.1,0.05,0.05,0.1),2);deltas<-c(0.29,0.44);t<-4;N=(t-1)*4;m<-26
  
rho2<-matrix(c(1,0.5,0.5,1),2)
  rho01<-matrix(c(0.01,0.005,0.005,0.01),2);deltas<-c(0.13,0.26);t<-4;N=(t-1)*5;m<-95
  rho01<-matrix(c(0.01,0.005,0.005,0.05),2);deltas<-c(0.17,0.32);t<-3;N=(t-1)*8;m<-92
  rho01<-matrix(c(0.01,0.005,0.005,0.1),2);deltas<-c(0.24,0.41);t<-5;N=(t-1)*5;m<-12
  
  rho01<-matrix(c(0.05,0.005,0.005,0.01),2);deltas<-c(0.43,0.3);t<-4;N=(t-1)*7;m<-10
  rho01<-matrix(c(0.05,0.025,0.025,0.05),2);deltas<-c(0.22,0.29);t<-5;N=(t-1)*4;m<-25
  rho01<-matrix(c(0.05,0.025,0.025,0.1),2);deltas<-c(0.15,0.5);t<-3;N=(t-1)*11;m<-83
  
  rho01<-matrix(c(0.1,0.005,0.005,0.01),2);deltas<-c(0.46,0.19);t<-3;N=(t-1)*13;m<-40
  rho01<-matrix(c(0.1,0.025,0.025,0.05),2);deltas<-c(0.5,0.1);t<-4;N=(t-1)*9;m<-88
  rho01<-matrix(c(0.1,0.05,0.05,0.1),2);deltas<-c(0.28,0.38);t<-5;N=(t-1)*6;m<-8

rho2<-matrix(c(1,0.8,0.8,1),2)
  rho01<-matrix(c(0.01,0.005,0.005,0.01),2);deltas<-c(0.2,0.11);t<-5;N=(t-1)*4;m<-89
  rho01<-matrix(c(0.01,0.005,0.005,0.05),2);deltas<-c(0.37,0.41);t<-3;N=(t-1)*15;m<-8
  rho01<-matrix(c(0.01,0.005,0.005,0.1),2);deltas<-c(0.35,0.4);t<-4;N=(t-1)*10;m<-5
  
  rho01<-matrix(c(0.05,0.005,0.005,0.01),2);deltas<-c(0.12,0.12);t<-5;N=(t-1)*5;m<-75
  rho01<-matrix(c(0.05,0.025,0.025,0.05),2);deltas<-c(0.5,0.5);t<-3;N=(t-1)*5;m<-22
  rho01<-matrix(c(0.05,0.025,0.025,0.1),2);deltas<-c(0.27,0.36);t<-4;N=(t-1)*6;m<-18
  
  rho01<-matrix(c(0.1,0.005,0.005,0.01),2);deltas<-c(0.41,0.15);t<-3;N=(t-1)*14;m<-60
  rho01<-matrix(c(0.1,0.025,0.025,0.05),2);deltas<-c(0.49,0.23);t<-4;N=(t-1)*8;m<-17
  rho01<-matrix(c(0.1,0.05,0.05,0.1),2);deltas<-c(0.32,0.23);t<-5;N=(t-1)*2;m<-56
  
calPower_IU(deltas,margins=c(0,0),vars=c(1,1),rho01,rho2,N,t,m,K=2,alpha=0.05)

#check
scenarios<-read.table("/Users/kdavis07/Dropbox/SW-CRT Methods Development/2_CoPrimary/RCode/Simulations/HH/Sim_Params.txt", header=TRUE, sep="")
power<-NULL
for(k in 1:27){
  scenario <- subset(scenarios, scenario == k)
  
  t <- scenario$t
  N <- scenario$N
  m <- scenario$m
  deltas<-c(scenario$delta1,scenario$delta2)
  rho01<-matrix(c(scenario$rho01.11,scenario$rho01.12,scenario$rho01.12,scenario$rho01.22),2)
  rho2<-matrix(c(1,scenario$rho2.12,scenario$rho2.12,1),2)
  
  pred<-calPower_IU(deltas,margins=c(0,0),vars=c(1,1),rho01,rho2,N,t,m,K=2,alpha=0.05)
  
  power.k <-cbind(pred$pred.power.t,pred$pred.power.z)
  
  power<-rbind(power,power.k)
}

# Naive
power<-NULL
for(k in 1:27){
  scenario <- subset(scenarios, scenario == k)
  
  t <- scenario$t
  N <- scenario$N
  m <- scenario$m
  deltas<-c(scenario$delta1,scenario$delta2)
  rho01<-matrix(c(scenario$rho01.11,0,0,scenario$rho01.22),2)
  rho2<-matrix(c(1,0,0,1),2)
  
  pred<-calPower_IU(deltas,margins=c(0,0),vars=c(1,1),rho01,rho2,N,t,m,K=2,alpha=0.05)
  
  power.k <-cbind(pred$pred.power.t,pred$pred.power.z)
  
  power<-rbind(power,power.k)
}








