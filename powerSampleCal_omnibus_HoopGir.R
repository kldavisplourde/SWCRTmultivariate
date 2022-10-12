library(mvtnorm)
########################################################################################################################################################
##Power/Sample Size Calculation based on the t-test for omnibus test.##
# INPUT
# deltas: (delta_1,...,delta_K), the vector of treatment effect for (1st,...,Kth) endpoints
# margins: (margin_1,...,margin_K), the vector of non-inferiority margins, when delta_1 = ... = delta_K = 0,
#          superiority tests are performed on all endpoints
# vars: (var_1,...,var_K), the vector of marginal variance for (1st,...,Kth) endpoints
# rho0: a K by K dimensional matrix for the correlation parameters (rho0^k) and (rho0^kk')
# For rho0: Endpoint-specific ICCs (within-period)
#           the diagonal elements correspond to rho0^k's 
#           the off-diagonal elements correspond to (rho0^kk')'s 
#           For example, rho0[1,1] corresponds to rho0^1, which is the within-period ICC for the first endpoint
#                        rho0[1,2] corresponds to rho0^12, which is the within-period correlation of outcomes between subjects on the 1st and 2nd endpoints    
# rho1: a K by K dimensional matrix for the correlation parameters (rho1^k) and (rho1^kk')
# For rho1: Endpoint-specific ICCs (between-period)
#           the diagonal elements correspond to rho1^k's 
#           the off-diagonal elements correspond to (rho1^kk')'s 
#           For example, rho1[1,1] corresponds to rho1^1, which is the between-period ICC for the first endpoint
#                        rho1[1,2] corresponds to rho1^12, which is the between-period correlation of outcomes between subjects on the 1st and 2nd endpoints    
# rho2: a K by K dimensional matrix for the correlation parameters (rho2^kk')
# For rho2: Intra-subject ICC
#           the diagonal elements are 1
#           the off-diagonal elements correspond to (rho2^kk')'s
#           For example, rho2[1,2] corresponds to rho2^12, which is the correlation of outcomes within same subject on the 1st and 2nd endpoints    
# N: number of clusters
# t: number of time periods
# m: cluster-period size
# K: number of endpoints
# alpha: upper bound of type I error rates over the whole null space
########################################################################################################################################################
####Function to Calculate Power Given Design Configurations based on the t test and normal (omnibus test)#######
####Critical values c_1,...,c_K are set to t_alpha, (1-alpha)th quantile of the t distribution with df = N-2K###
calPower_omnibus <- function(deltas,margins,vars,rho0,rho1,rho2,N,t,m,K,alpha)
{
  # Common treatment effects?
  if(length(deltas)==1){common.treat.eff <- 1} else {common.treat.eff <- 0}
  
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
  constrRiE <- function(rho0,rho2,K,vars)
  { rho0k <- diag(rho0)
  SigmaE_Matrix <- diag((1-rho0k)*vars)
  for(row in 1:K )
  {
    for(col in 1:K)
    {
      if(row != col){
        SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]-rho0[row,col])
      }
    }
  }
  return(SigmaE_Matrix)
  }
  
  #####function to construct covariance matrix Sigma_phi for Y_i########
  constrRiP <- function(rho1,K,vars)
  { rho0k <- diag(rho1)
  SigmaP_Matrix <- diag(rho0k*vars)
  for(row in 1:K )
  {
    for(col in 1:K)
    {
      if(row != col){
        SigmaP_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho1[row,col]
      }
    }
  }
  return(SigmaP_Matrix)
  }
  
  #####function to construct covariance matrix Sigma_psi for Y_i########
  constrRiPs <- function(rho1,rho0,K,vars)
  { rho0k <- diag(rho0) - diag(rho1)
  SigmaPs_Matrix <- diag(rho0k*vars)
  for(row in 1:K )
  {
    for(col in 1:K)
    {
      if(row != col){
        SigmaPs_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho0[row,col]-rho1[row,col])
      }
    }
  }
  return(SigmaPs_Matrix)
  }
  
  ####Define function to calculate covariance between deltas#####
  calCovbetas <- function(vars,rho1,rho0,rho2){
    sigmaE <- constrRiE(rho0,rho2,K,vars)
    sigmaP <- constrRiP(rho1,K,vars)
    sigmaPs <- constrRiPs(rho1,rho0,K,vars)
    #Stepped wedge
    if(common.treat.eff == 1){
      sde <- matrix(c(sqrt(sigmaE[1,1]),sqrt(sigmaE[2,2])),2,1)
      covMatrix <- N*t*solve((N*t*U-t*W+U^2-N*V)*t(sde)%*%solve(sigmaPs +(1/m)*sigmaE)%*%sde-(U^2-N*V)*t(sde)%*%solve(t*sigmaP+sigmaPs+(1/m)*sigmaE)%*%sde)
    } else {
      #covMatrix <- ((N*t)/(N*t*U-t*W+U^2-N*V))*solve(solve(sigmaPs+(1/m)*sigmaE)+((N*V-U^2)/(N*t*U-t*W+U^2-N*V))*solve(t*sigmaP+sigmaPs+(1/m)*sigmaE))
      covMatrix <- N*t*solve((N*t*U-t*W+U^2-N*V)*solve(sigmaPs +(1/m)*sigmaE)-(U^2-N*V)*solve(t*sigmaP+sigmaPs+(1/m)*sigmaE))
    }
    
    return(covMatrix)  
  }
  
  ####Define function to calculate correlation between test statistics #####
  calCorWks <-  function(vars,rho1,rho0,rho2)
  {
    top <- calCovbetas(vars,rho1,rho0,rho2)
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
  sigmaks.sq <- diag(calCovbetas(vars,rho1,rho0,rho2))
  meanVector <- (deltas-margins)/sqrt(sigmaks.sq)
  
  if(common.treat.eff == 1){
    criticalValue.t <- qt(p=(1-alpha/2), df=(N-K-1))
    pred.power <- 1-2*pt(criticalValue.t, df=(N-K-1), meanVector)[1]
  } else {
    omega <- calCovbetas(vars,rho1,rho0,rho2)
    #Variance Inflation of multivariate outcomes
    tau <- N*t(deltas) %*% solve(omega) %*% deltas
    Fscore <- qf(1-alpha, df1=K, df2=(N-2*K), ncp=0, lower.tail = TRUE, log.p = FALSE)
    pred.power <-1- pf(Fscore, df1=K, df2=(N-2*K), ncp=tau, lower.tail = TRUE, log.p = FALSE)
  }
  
  param <- list(vard=c(sigmaks.sq),pred.power=pred.power)
  return(param)
}