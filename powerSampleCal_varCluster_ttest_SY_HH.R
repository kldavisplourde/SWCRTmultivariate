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
  ####Define function to calculate covariance between betas#####
  calCovbetas <- function(vars,rho01,rho2){
    sigmaE <- constrRiE(rho01,rho2,K,vars)
    sigmaP <- constrRiP(rho01,K,vars)
    #tmp <- solve(diag(1,K)-cv^2*(m*sigmaP %*% solve(sigmaE + m*sigmaP) %*% sigmaE %*% solve(sigmaE + m*sigmaP) ))
    #covMatrix <- 1/(m*sigmaz.square)*(sigmaE+m*sigmaP)%*%tmp
    #covMatrix <- (covMatrix +t(covMatrix))/2  # symmerize the off-diagonal
    covMatrix <- (N/(m*(N*U-W)))*(sigmaE-((U^2-N*V)/(U^2+N*t*U-t*W-N*V))*solve(solve(sigmaE)+((N*U-W)/(m*(U^2+N*t*U-t*W-N*V)))*solve(sigmaP)))
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
  meanVector <- sqrt(N)*(deltas-margins)/sqrt(sigmaks.sq)
  wCor <- calCorWks(vars,rho01,rho2)
  criticalValue.t <- qt(p=(1-alpha), df=(N-2*K))
  criticalValue.n <- qnorm(p=(1-alpha))
  pred.power.t <- pmvt(lower = rep(criticalValue.t,K),upper=rep(Inf,K),df = (N-2*K), sigma = wCor,delta=meanVector)[1]
  pred.power.n <- pmvnorm(lower = rep(criticalValue.n,K),upper=rep(Inf,K), sigma = wCor,mean=meanVector)[1]
  
  param <- list(pred.power.t=pred.power.t,pred.power.n=pred.power.n)
  return(param)
}

deltas<-c(0.08,0.13); margins<-c(0,0); vars<-c(1.1,1.05); rho01<-matrix(c(0.02,0.01,0.01,0.015),2); rho2<-matrix(c(1,0.05,0.05,1),2); t<-4; N=(t-1)*5; m<-15;K<-2; alpha<-0.05
deltas<-c(0.8,0.045); margins<-c(0,0); vars<-c(1.5,1.2); rho01<-matrix(c(0.2,0.1,0.1,0.15),2); rho2<-matrix(c(1,0.3,0.3,1),2); t<-5; N=(t-1)*7; m<-12;K<-2; alpha<-0.05
deltas<-c(0.25,0.6); margins<-c(0,0); vars<-c(0.8,1.5); rho01<-matrix(c(0.04,0.02,0.02,0.05),2); rho2<-matrix(c(1,0.09,0.09,1),2); t<-3; N=(t-1)*4; m<-10;K<-2; alpha<-0.05
deltas<-c(1.2,0.2); margins<-c(0,0); vars<-c(2,1.5); rho01<-matrix(c(0.06,0.02,0.02,0.04),2); rho2<-matrix(c(1,0.1,0.1,1),2); t<-3; N=(t-1)*8; m<-5;K<-2; alpha<-0.05
deltas<-c(0.22,0.22); margins<-c(0,0); vars<-c(1.1,0.7); rho01<-matrix(c(0.1,0.05,0.05,0.15),2); rho2<-matrix(c(1,0.2,0.2,1),2); t<-4; N=(t-1)*4; m<-4;K<-2; alpha<-0.05

calPower_IU(deltas,margins,vars,rho01,rho2,N,t,m,K,alpha)
