library(mvtnorm)
########################################################################################################################################################
##Power/Sample Size Calculation based on the t test##
# INPUT
# deltas: (delta_1,...,delta_K), the vector of treatment effect for (1st,...,Kth) endpoints
# margins: (margin_1,...,margin_K), the vector of non-inferiority margins, when delta_1 = ... = delta_K = 0,
#          superiority tests are performed on all endpoints
# vars: (var_1,...,var_K), the vector of marginal variance for (1st,...,Kth) endpoints
# rho01: a K by K dimensional matrix for the correlation parameters (rho0^k) and (rho1^kk')
# For rho01: cluster random effect
#           the diagonal elements correspond to rho0^k's 
#           the off-diagonal elements correspond to (rho1^kk')'s 
#           For example, rho01[1,1] corresponds to rho0^1, which is the ICC for the first endpoint
#                        rho01[1,2] corresponds to rho1^12, which is the correlation of outcomes between subjects on the 1st and 2nd endpoints    
# rho02: a K by K dimensional matrix for the correlation parameters (rho0^k) and (rho1^kk')
# For rho02: cluster-period random effect
#           the diagonal elements correspond to rho0^k's 
#           the off-diagonal elements correspond to (rho1^kk')'s 
#           For example, rho02[1,1] corresponds to rho0^1, which is the ICC for the first endpoint
#                        rho02[1,2] corresponds to rho1^2, which is the correlation of outcomes between subjects on the 1st and 2nd endpoints    
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
calPower_IU <- function(deltas,margins,vars,rho01,rho02,rho2,N,t,m,K,alpha)
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
  constrRiE <- function(rho02,rho2,K,vars)
  { rho0k <- diag(rho02)
    SigmaE_Matrix <- diag((1-rho0k)*vars)
    for(row in 1:K )
    {
      for(col in 1:K)
      {
        if(row != col){
          SigmaE_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]-rho02[row,col])
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
  
  #####function to construct covariance matrix Sigma_psi for Y_i########
  constrRiPs <- function(rho01,rho02,K,vars)
  { rho0k <- diag(rho02) - diag(rho01)
  SigmaPs_Matrix <- diag(rho0k*vars)
  for(row in 1:K )
  {
    for(col in 1:K)
    {
      if(row != col){
        SigmaPs_Matrix[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho02[row,col]-rho01[row,col])
      }
    }
  }
  return(SigmaPs_Matrix)
  }
  
  ####Define function to calculate covariance between deltas#####
  calCovbetas <- function(vars,rho01,rho02,rho2){
    sigmaE <- constrRiE(rho02,rho2,K,vars)
    sigmaP <- constrRiP(rho01,K,vars)
    sigmaPs <- constrRiPs(rho01,rho02,K,vars)
  #Stepped wedge
    #covMatrix <- ((N*t)/(N*t*U-t*W+U^2-N*V))*solve(solve(sigmaPs+(1/m)*sigmaE)+((N*V-U^2)/(N*t*U-t*W+U^2-N*V))*solve(t*sigmaP+sigmaPs+(1/m)*sigmaE))
    covMatrix <- N*t*solve((N*t*U-t*W+U^2-N*V)*solve(sigmaPs +(1/m)*sigmaE)-(U^2-N*V)*solve(t*sigmaP+sigmaPs+(1/m)*sigmaE))
    return(covMatrix)  
  }
  
  ####Define function to calculate correlation between test statistics #####
  calCorWks <-  function(vars,rho01,rho02,rho2)
  {
    top <- calCovbetas(vars,rho01,rho02,rho2)
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
  sigmaks.sq <- diag(calCovbetas(vars,rho01,rho02,rho2))
  meanVector <- (deltas-margins)/sqrt(sigmaks.sq)
  wCor <- calCorWks(vars,rho01,rho02,rho2)
  criticalValue.t <- qt(p=(1-alpha), df=(N-2*K))
  criticalValue.n <- qnorm(p=(1-alpha))
  pred.power.t <- pmvt(lower = rep(criticalValue.t,K),upper=rep(Inf,K),df = (N-2*K), sigma = wCor,delta=meanVector)[1]
  pred.power.n <- pmvnorm(lower = rep(criticalValue.n,K),upper=rep(Inf,K), sigma = wCor,mean=meanVector)[1]
  
  param <- list(vard=c(sigmaks.sq),pred.power.t=pred.power.t,pred.power.z=pred.power.n)
  return(param)
}

#Sim Study Development
#sd<-2
#rho2<-matrix(c(1,0.2,0.2,1),2)
  #rho02<-matrix(c(0.02,0.01,0.01,0.02),2);rho01<-matrix(c(0.01,0.005,0.005,0.01),2);deltas<-c(0.86*sd,0.71*sd);t<-3;N=(t-1)*4;m<-13
  #rho02<-matrix(c(0.02,0.01,0.01,0.1),2);rho01<-matrix(c(0.01,0.005,0.005,0.05),2);deltas<-c(0.6*sd,0.78*sd);t<-5;N=(t-1)*2;m<-8
  #rho02<-matrix(c(0.02,0.01,0.01,0.2),2);rho01<-matrix(c(0.01,0.005,0.005,0.1),2);deltas<-c(0.45*sd,0.92*sd);t<-4;N=(t-1)*3;m<-15
  
  #rho02<-matrix(c(0.1,0.01,0.01,0.02),2);rho01<-matrix(c(0.05,0.005,0.005,0.01),2);deltas<-c(0.82*sd,0.65*sd);t<-5;N=(t-1)*2;m<-7
  #rho02<-matrix(c(0.1,0.05,0.05,0.1),2);rho01<-matrix(c(0.05,0.025,0.025,0.05),2);deltas<-c(0.49*sd,0.98*sd);t<-4;N=(t-1)*4;m<-15
  #rho02<-matrix(c(0.1,0.05,0.05,0.2),2);rho01<-matrix(c(0.05,0.025,0.025,0.1),2);deltas<-c(0.59*sd,0.99*sd);t<-3;N=(t-1)*6;m<-20
  
  #rho02<-matrix(c(0.2,0.01,0.01,0.02),2);rho01<-matrix(c(0.1,0.005,0.005,0.01),2);deltas<-c(0.47*sd,0.22*sd);t<-5;N=(t-1)*5;m<-18
  #rho02<-matrix(c(0.2,0.05,0.05,0.1),2);rho01<-matrix(c(0.1,0.025,0.025,0.05),2);deltas<-c(0.92*sd,0.92*sd);t<-3;N=(t-1)*5;m<-12
  #rho02<-matrix(c(0.2,0.1,0.1,0.2),2);rho01<-matrix(c(0.1,0.05,0.05,0.1),2);deltas<-c(0.54*sd,0.81*sd);t<-4;N=(t-1)*4;m<-25

#rho2<-matrix(c(1,0.5,0.5,1),2)
  #rho02<-matrix(c(0.02,0.01,0.01,0.02),2);rho01<-matrix(c(0.01,0.005,0.005,0.01),2);deltas<-c(0.62*sd,0.53*sd);t<-4;N=(t-1)*5;m<-5
  #rho02<-matrix(c(0.02,0.01,0.01,0.1),2);rho01<-matrix(c(0.01,0.005,0.005,0.05),2);deltas<-c(0.34*sd,0.88*sd);t<-3;N=(t-1)*8;m<-22
  #rho02<-matrix(c(0.02,0.01,0.01,0.2),2);rho01<-matrix(c(0.01,0.005,0.005,0.1),2);deltas<-c(0.24*sd,0.83*sd);t<-5;N=(t-1)*5;m<-13
  
  #rho02<-matrix(c(0.1,0.01,0.01,0.02),2);rho01<-matrix(c(0.05,0.005,0.005,0.01),2);deltas<-c(0.38*sd,0.55*sd);t<-4;N=(t-1)*7;m<-10
  #rho02<-matrix(c(0.1,0.05,0.05,0.1),2);rho01<-matrix(c(0.05,0.025,0.025,0.05),2);deltas<-c(0.3*sd,0.68*sd);t<-5;N=(t-1)*4;m<-25
  #rho02<-matrix(c(0.1,0.05,0.05,0.2),2);rho01<-matrix(c(0.05,0.025,0.025,0.1),2);deltas<-c(0.62*sd,0.62*sd);t<-3;N=(t-1)*11;m<-8

  #rho02<-matrix(c(0.2,0.01,0.01,0.02),2);rho01<-matrix(c(0.1,0.005,0.005,0.01),2);deltas<-c(0.84*sd,0.29*sd);t<-3;N=(t-1)*13;m<-18
  #rho02<-matrix(c(0.2,0.05,0.05,0.1),2);rho01<-matrix(c(0.1,0.025,0.025,0.05),2);deltas<-c(0.4*sd,0.4*sd);t<-4;N=(t-1)*9;m<-11
  #rho02<-matrix(c(0.2,0.1,0.1,0.2),2);rho01<-matrix(c(0.1,0.05,0.05,0.1),2);deltas<-c(0.32*sd,0.84*sd);t<-5;N=(t-1)*6;m<-24
  
#rho2<-matrix(c(1,0.8,0.8,1),2)
  #rho02<-matrix(c(0.02,0.01,0.01,0.02),2);rho01<-matrix(c(0.01,0.005,0.005,0.01),2);deltas<-c(0.24*sd,0.55*sd);t<-5;N=(t-1)*4;m<-20
  #rho02<-matrix(c(0.02,0.01,0.01,0.1),2);rho01<-matrix(c(0.01,0.005,0.005,0.05),2);deltas<-c(0.29*sd,0.57*sd);t<-3;N=(t-1)*15;m<-14
  #rho02<-matrix(c(0.02,0.01,0.01,0.2),2);rho01<-matrix(c(0.01,0.005,0.005,0.1),2);deltas<-c(0.2*sd,0.84*sd);t<-4;N=(t-1)*10;m<-17

  #rho02<-matrix(c(0.1,0.01,0.01,0.02),2);rho01<-matrix(c(0.05,0.005,0.005,0.01),2);deltas<-c(0.31*sd,0.62*sd);t<-5;N=(t-1)*5;m<-13
  #rho02<-matrix(c(0.1,0.05,0.05,0.1),2);rho01<-matrix(c(0.05,0.025,0.025,0.05),2);deltas<-c(0.63*sd,0.92*sd);t<-3;N=(t-1)*5;m<-22
  #rho02<-matrix(c(0.1,0.05,0.05,0.2),2);rho01<-matrix(c(0.05,0.025,0.025,0.1),2);deltas<-c(0.45*sd,0.45*sd);t<-4;N=(t-1)*6;m<-18

  #rho02<-matrix(c(0.2,0.01,0.01,0.02),2);rho01<-matrix(c(0.1,0.005,0.005,0.01),2);deltas<-c(0.99*sd,0.25*sd);t<-3;N=(t-1)*14;m<-25
  #rho02<-matrix(c(0.2,0.05,0.05,0.1),2);rho01<-matrix(c(0.1,0.025,0.025,0.05),2);deltas<-c(0.63*sd,0.31*sd);t<-4;N=(t-1)*8;m<-17
  #rho02<-matrix(c(0.2,0.1,0.1,0.2),2);rho01<-matrix(c(0.1,0.05,0.05,0.1),2);deltas<-c(0.82*sd,0.82*sd);t<-5;N=(t-1)*2;m<-10

#calPower_IU(deltas,margins=c(0,0),vars=c(4,4),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)

#check
scenarios<-read.table("/Users/kdavis07/Dropbox/SW-CRT Methods Development/2_CoPrimary/RCode/Simulations/HoopGir/Sim_Params.txt", header=TRUE, sep="")
power<-NULL
for(k in 1:27){
  scenario <- subset(scenarios, scenario == k)
  
  t <- scenario$t
  N <- scenario$N
  m <- scenario$m
  deltas<-c(scenario$delta1,scenario$delta2)
  rho01<-matrix(c(scenario$rho01.11,scenario$rho01.12,scenario$rho01.12,scenario$rho01.22),2)
  rho02<-matrix(c(scenario$rho02.11,scenario$rho02.12,scenario$rho02.12,scenario$rho02.22),2)
  rho2<-matrix(c(1,scenario$rho2.12,scenario$rho2.12,1),2)
  
  pred<-calPower_IU(deltas,margins=c(0,0),vars=c(4,4),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
  
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
  rho02<-matrix(c(scenario$rho02.11,0,0,scenario$rho02.22),2)
  rho2<-matrix(c(1,0,0,1),2)
  
  pred<-calPower_IU(deltas,margins=c(0,0),vars=c(4,4),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
  
  power.k <-cbind(pred$pred.power.t,pred$pred.power.z)
  
  power<-rbind(power,power.k)
}








