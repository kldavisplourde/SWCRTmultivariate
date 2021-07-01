#setwd("C:/Users/lifan/Dropbox/CRT_Coprimary/Data/")
#data=read.csv("short2.csv",header=TRUE)

#library(nlme)
library(lme4)
library(mvtnorm)
library(numDeriv)

# function to perform EM estimation with K=2 outcomes
EM.estim <- function(data, fm1,fm2, cluster,cluster.period, maxiter=500,epsilon=1e-4
                     , verbose=FALSE){
  # fit mixed model to initialize parameters
  #fm1 <- lme(formula1, random = ~ 1|cluster, data=data,control=lmeControl(returnObject=TRUE))
  #fm2 <- lme(formula2, random = ~ 1|cluster, data=data,control=lmeControl(returnObject=TRUE))
  K <- 2
  zeta <- as.numeric(c(fixef(fm1), fixef(fm2)))
  beta1 = as.numeric(fixef(fm1))
  beta2 = as.numeric(fixef(fm2))
  if (length(beta1) != length(beta2))
    stop("\nnumber of covariates do not match between endpoints.")
  nvar<-length(beta1)
  TermsX1 <- terms(fm1)
  TermsX2 <- terms(fm2)
  mfX1 <- model.frame(TermsX1, data = data)[,-1]
  mfX2 <- model.frame(TermsX2, data = data)[,-1]
  if (identical(mfX1,mfX2) == FALSE)
    stop("\ncovariates do not match between endpoints.")
  ##
  
  # vector of cluster sizes and cluster-period sizes
  m <- table(data[, paste(cluster)])
  mp <- table(data[, paste(cluster.period)])
  
  vc1<-as.data.frame(VarCorr(fm1))
  vc2<-as.data.frame(VarCorr(fm2))
  s2phi1 <- vc1[vc1$grp==paste(cluster),4]
  s2phi2 <- vc2[vc2$grp==paste(cluster),4]
  SigmaPhi <- diag(c(s2phi1, s2phi2))
  InvS2Phi <- solve(SigmaPhi)
  
  s2psi1 <- vc1[vc1$grp==paste(cluster.period),4]
  s2psi2 <- vc2[vc2$grp==paste(cluster.period),4]
  SigmaPsi <- diag(c(s2psi1, s2psi2))
  InvS2Psi <- solve(SigmaPsi)
  
  s2e1 <- vc1[vc1$grp=="Residual",4]
  s2e2 <- vc2[vc2$grp=="Residual",4]
  SigmaE <- diag(c(s2e1, s2e2))
  InvS2E <- solve(SigmaE)
  
  Y <- as.matrix(cbind(model.frame(TermsX1, data = data)[,1],model.frame(TermsX2, data = data)[,1]))
  ID <- data[, paste(cluster)]
  ID.period <- data[, paste(cluster.period)]
  n <- length(unique(ID))
  np <- length(unique(data[, paste(cluster.period)]))
  nperiods <- table(unique(cbind(data[, paste(cluster)],data[, paste(cluster.period)]))[,1])
  cID.period <- unique(cbind(data[, paste(cluster)],data[, paste(cluster.period)]))[,1]
  X <- as.matrix(cbind(1, mfX1)) # design matrix
  
  ESSphi1 <- matrix(0,n,K)
  ESSphi2 <- array(0,c(K,K,n))
  
  ESSpsi1 <- matrix(0,np,K)
  ESSpsi2 <- array(0,c(K,K,np))
  
  phipsi2 <- array(0,c(K,K,np))
  ESSphi2vpsi <- array(0,c(K,K,n))
  Vjpsi <- array(0,c(K,K,np))
  ESSsumpsi2 <- array(0,c(K,K,n))
  ESSpsi1n2 <- matrix(0,np,K)
  
  
  #maxiter=500
  #epsilon=1e-4
  delta = 2*epsilon
  max_modi = 20
  
  converge = 0
  
  # log likelihood
  
  loglik = function(theta){
    beta1 = theta[1:nvar]
    beta2 = theta[(nvar+1):(2*nvar)]
    sphi11 = theta[(2*nvar+1)]
    sphi12 = theta[(2*nvar+2)]
    sphi22 = theta[(2*nvar+3)]
    se11 = theta[(2*nvar+4)]
    se12 = theta[(2*nvar+5)]
    se22 = theta[(2*nvar+6)]
    spsi11 = theta[(2*nvar+7)]
    spsi12 = theta[(2*nvar+8)]
    spsi22 = theta[(2*nvar+9)]
    SigmaPhi = matrix(c(sphi11,sphi12,sphi12,sphi22),2,2)
    SigmaPsi = matrix(c(spsi11,spsi12,spsi12,spsi22),2,2)
    SigmaE = matrix(c(se11,se12,se12,se22),2,2)
    InvS2Phi <- solve(SigmaPhi)
    InvS2Psi <- solve(SigmaPsi)
    InvS2E <- solve(SigmaE)
    
    temp <- 0
    for(j in 1:n){
      N <- m[j]/nperiods[j]
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
      obs = c(t(residj))
      tm1 <- nperiods[j]*(N-1)*log(det(SigmaE))+(nperiods[j]-1)*log(det(SigmaE+N*SigmaPsi))+log(det(SigmaE+N*(SigmaPsi+nperiods[j]*SigmaPhi)))
      InvSS2 <- solve(SigmaE+N*SigmaPsi)-InvS2E
      InvSS22 <- solve(SigmaE+N*(SigmaPsi+nperiods[j]*SigmaPhi))-solve(SigmaE+N*SigmaPsi)
      Invj <- kronecker(diag(nrow=nperiods[j]),kronecker(diag(nrow=N),InvS2E) + 
        kronecker(matrix(1,N,N),InvSS2)/N) +
        kronecker(matrix(1,nperiods[j],nperiods[j]),kronecker(matrix(1,N,N),InvSS22)/(nperiods[j]*N))
      tm2 <- c(t(obs) %*% Invj %*% obs)
      temp <- temp-(tm1+tm2)/2
    }
    temp
  }
  thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]),c(SigmaPsi[!lower.tri(SigmaPsi)]))
  LLold <- loglik(thetah)
  
  
  niter=1
  while((niter <= maxiter) & (abs(delta) > epsilon)){
    
    # Expectation step
    for(j in 1:np){
      Yj <- Y[ID.period == j,,drop=FALSE]
      Xj <- X[ID.period == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
      Vj <- solve(InvS2Psi + mp[j]*InvS2E)
      Muj <- as.numeric(Vj %*% InvS2E %*% colSums(residj))
      Nujj <- Vj + tcrossprod(Muj)
      ESSpsi1[j,] <- Muj
      ESSpsi2[,,j] <- Nujj
      
      ESSpsi1n2[j,] <- (mp[j]*mp[j])*Muj
      Vjpsi[,,j] <- (mp[j]*mp[j])*Vj
      phipsi2[,,j] <- (mp[j]*mp[j])*(Vj %*% InvS2E)
    }
    
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
      Vj <- solve(InvS2Phi + m[j]*InvS2E)
      Muj <- as.numeric(Vj %*% InvS2E %*% colSums(residj))
      Nujj <- Vj + tcrossprod(Muj)
      ESSphi1[j,] <- Muj
      ESSphi2[,,j] <- Nujj
      
      ind.f <- min(which(cID.period == j))
      ind.l <- max(which(cID.period == j))
      phipsi2j <- phipsi2[,,ind.f:ind.l]
      ESSphi2vpsi[,,j] <- Nujj%*%rowSums(phipsi2j,dims=2)
      
      Vpsij <- Vjpsi[,,ind.f:ind.l]
      ESSsumpsi2[,,j] <- (Vj %*% InvS2E)%*%(rowSums(Vpsij,dims=2) + tcrossprod(colSums(ESSpsi1n2[ind.f:ind.l,])))
    }
    
    # Maximization step - phi & psi
    SigmaPhi <- apply(ESSphi2,1:2, sum)/n
    InvS2Phi <- solve(SigmaPhi)
    
    SigmaPsi <- apply(ESSpsi2,1:2, sum)/np
    InvS2Psi <- solve(SigmaPsi)
    
    # Maximization step - zeta
    # Simplify the expression analytically, and obtain simple expression!
    XXt <- crossprod(X)
    Vzeta <-solve(kronecker(InvS2E, XXt))
    rzeta1 <- t(X)%*%(Y[,1]-ESSphi1[ID,1]-ESSpsi1[ID.period,1])
    rzeta2 <- t(X)%*%(Y[,2]-ESSphi1[ID,2]-ESSpsi1[ID.period,2])
    zeta <- Vzeta %*% rbind(InvS2E[1,1]*rzeta1 + InvS2E[1,2]*rzeta2,
                            InvS2E[2,1]*rzeta1 + InvS2E[2,2]*rzeta2)
    zeta <- c(zeta)
    beta1 = zeta[1:nvar]
    beta2 = zeta[(nvar+1):(2*nvar)]
    
    # Maximization step - epsilon
    re <- Y - cbind(X%*%beta1, X%*%beta2)
    rss <- crossprod(re) + rowSums(sweep(ESSphi2,3,m,FUN="*"),dims=2) + rowSums(sweep(ESSpsi2,3,mp,FUN="*"),dims=2) -
      crossprod(ESSphi1,rowsum(re,ID)) - crossprod(rowsum(re,ID),ESSphi1) -
      crossprod(ESSpsi1,rowsum(re,ID.period)) - crossprod(rowsum(re,ID.period),ESSpsi1) +
      crossprod(ESSphi1,rowsum(sweep(ESSpsi1,1,mp,FUN="*"),cID.period)) + crossprod(rowsum(sweep(ESSpsi1,1,mp,FUN="*"),cID.period),ESSphi1) +
      rowSums(ESSphi2vpsi,dims=2) +rowSums(ESSsumpsi2,dims=2)
    SigmaE <- rss/sum(m)
    InvS2E <- solve(SigmaE)
    
    # whether the algorithm converges
    thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]),c(SigmaPsi[!lower.tri(SigmaPsi)]))
    LLnew <- loglik(thetah)
    delta <- abs(LLnew - LLold)
    LLold <- LLnew
    converge = (abs(delta)<=epsilon)
    niter <- niter + 1
    #if(verbose) cat(paste('iter=',niter),'\t',paste('param.error=',epsilon),'\t',paste('loglik=',LLnew),'\n');
    if(verbose) cat(paste('iter=',niter),'\t',paste('param.error=',epsilon),'\t',paste('loglik=',LLnew),'\t',paste('SigmaE=',c(SigmaE[!lower.tri(SigmaE)])),'\n')
    
    #print(niter)
    #print(zeta)
    #print(SigmaPhi)
    #print(SigmaE)
    #print(LLnew)
  }
  param <- list(theta=list(zeta=zeta,SigmaE=SigmaE,SigmaPhi=SigmaPhi,SigmaPsi=SigmaPsi),loglik=LLnew,eps=epsilon,iter=niter)
  return(param) 
}







