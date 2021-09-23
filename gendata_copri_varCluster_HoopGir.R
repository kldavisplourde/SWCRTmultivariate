# Original code from Fan Li, modified by Kendra Plourde

###############################################
### function for generating the SW-CRT data ###
# eff: (delta_1,...,delta_K), the vector of treatment effect for (1st,...,Kth) endpoints
# time.eff: (beta_11,...,beta_T1,...,beta_1K,...,beta_TK), the vector of time effects for (1st,...,Kth) endpoints
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
# n: number of clusters
# m: cluster size
# K: number of endpoints (K=2)
###############################################
datagen_cont <- function(n, m, K, cv, rho01, rho02, rho2, vars, eff, time.eff ){
  
  library(MASS)
  
  #####function to construct covariance matrix Sigma_E for Y_i########
  rho0k <- diag(rho02)
  sigmae <- diag((1-rho0k)*vars)
  for(row in 1:K )
  {
    for(col in 1:K)
    {
      if(row != col){
        sigmae[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho2[row,col]-rho02[row,col])
      }
    }
  }
  
  #####function to construct covariance matrix Sigma_phi for Y_i########
  rho0k <- diag(rho01)
  sigmac <- diag(rho0k*vars)
  for(row in 1:K )
  {
    for(col in 1:K)
    {
      if(row != col){
        sigmac[row,col] <- sqrt(vars[row])*sqrt(vars[col])*rho01[row,col]
      }
    }
  }
  
  #####function to construct covariance matrix Sigma_psi for Y_i########
  rho0k <- diag(rho02) - diag(rho01)
  sigmacp <- diag(rho0k*vars)
  for(row in 1:K )
  {
    for(col in 1:K)
    {
      if(row != col){
        sigmacp[row,col] <- sqrt(vars[row])*sqrt(vars[col])*(rho02[row,col]-rho01[row,col])
      }
    }
  }
  
  DATA<-NULL
  temp<-NULL 
  nperiods<-length(time.eff)/2 + 1
  if (cv>0){
  #parameters in gamma distribution
  kappa <- 1/cv^2
  theta <- m*cv^2
  
  # generate m from gamma
  mj <- round(rgamma(n, shape=kappa,  scale = theta))

 
  }else{
    #total sample size
    mj <- rep(m, n)
  }
  #total sample size
  #n=sum(mj)
  n.sub=sum(mj)*nperiods                               # KP: Assuming cluster size does not change over time. n = I*N*T in notes
  # number of cluster-periods
  mjp <- rep(mj, each=nperiods)                        # KP: Assuming cluster size does not change over time.
  #cluster <- rep(1:(2*Narm),  mj);
  cluster <- rep(1:n,  mj*nperiods)
  cluster.period <- rep(1:(n*nperiods),  mjp)
  
  
  #generate patient id
  patid <- seq(1:n)
  temp <- cbind(patid, cluster, cluster.period)
  temp <- as.data.frame(temp)
  
  #randomize treatment assignment
  #T <- rep(0, 2*Narm)
  #t <- sample(x=1:(2*Narm), size=Narm)
  #T[t] <- 1
  trtSeq<-matrix(0,nperiods-1,nperiods)
  trtSeq[upper.tri(trtSeq)]<-1
  trtSeq<-trtSeq[rep(1:nrow(trtSeq), rep(n/(nperiods-1),(nperiods-1))),]
  if(sum(rep(n/(nperiods-1),(nperiods-1))) != n){stop("invalid randomization scheme")}
  
  # time: need to repeat each element in c(0,1,...,T) mj times and concatenate all clusters
  for(i in 1:length(mj)) {if(i==1){time<-rep(0:(nperiods-1), each=mj[1])}else{if(i>1){time<-c(time,rep(0:(nperiods-1), each=mj[i]))}}}
  temp2 <- cbind(temp, time)
  #arm=0 if in the control group, arm=1 if in the intervention group
  #arm <- rep(T,  mj)
  temp2$arm<-0
  for(i in 1:nrow(temp2)) {
    temp2$arm[i]<-trtSeq[temp2$cluster[i],(temp2$time[i]+1)]
  }
  # create indicator variables for time
  for(i in 1:(nperiods-1)) {temp2[, paste0("time.", i)]<-ifelse(temp2$time==i,1,0)}
  
  DATA <- temp2[rep(1:nrow(temp2), each = K), ]
  
  #individual random effect
  # Parameters for bivariate normal distribution
  mu <- rep(0,K) # Mean

  epsilon <- mvrnorm(n*nperiods*m, mu = mu, Sigma = sigmae ) # from MASS package
  
  
  #cluster random effect
  g <- mvrnorm(n, mu = mu, Sigma = sigmac )
  #replicate row
  gamma <- g[rep(1:nrow(g),  mj*nperiods), ]

  #cluster-period random effect
  gp <- mvrnorm(n*nperiods, mu = mu, Sigma = sigmacp )
  #replicate row
  gammap <- gp[rep(1:nrow(gp),  mjp), ]
  
  # time effects
  time.effect<-rep(0,K)
  for(i in 1:(nperiods-1)){
    time.effect<-time.effect+matrix(temp2[, paste0("time.", i)]) %*% matrix(time.eff[c(i,(i+nperiods-1))],nrow=1) 
  }
  Y <- time.effect + matrix(temp2$arm) %*% matrix(eff,nrow=1) +gamma+gammap+epsilon
  colnames(Y) <- paste("out",seq(1,K),sep="")
  Yc <- as.vector(t(Y))
  label = rep(1:K,n)
  # design matrix
  full <- cbind(DATA, Yc, label)
  full <- data.frame(full)
  short <- cbind(temp2,Y )
  return(list(long=full, short=short))
  }
  
 