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
  # Z matrix
  bar.f <- findbars(formula(fm1)) # Identify random effect terms (find |)
  mf <- model.frame(subbars(fm1),data=data) # Replaces | with +
  Z <- t(mkReTrms(bar.f,mf,reorder.terms=FALSE)$Zt) # phi=(b_{11},b_{12},...,b_{I1},b_{I2},s_{111},s_{112},...,s_{1T1},s_{1T2},s_{I11},...,s_{IT2})'
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
  n <- length(unique(ID)) # number of clusters
  np <- length(unique(data[, paste(cluster.period)])) # number of cluster-periods
  nperiods <- table(unique(cbind(data[, paste(cluster)],data[, paste(cluster.period)]))[,1]) # number of periods in each cluster
  t <- max(nperiods)
  cID.period <- unique(cbind(data[, paste(cluster)],data[, paste(cluster.period)]))[,1]
  X <- as.matrix(cbind(1, mfX1)) # design matrix
  
  # Function for generating inverse of block diagonal covariance matrix for the set of all random effects
  bdiag_m <- function(lmat) {
    ## Copyright (C) 2016 Martin Maechler, ETH Zurich
    if(!length(lmat)) return(new("dgCMatrix"))
    stopifnot(is.list(lmat), is.matrix(lmat[[1]]),
              (k <- (d <- dim(lmat[[1]]))[1]) == d[2], # k x k
              all(vapply(lmat, dim, integer(2)) == k)) # all of them
    N <- length(lmat)
    if(N * k > .Machine$integer.max)
      stop("resulting matrix too large; would be  M x M, with M=", N*k)
    M <- as.integer(N * k)
    ## result: an   M x M  matrix
    new("dgCMatrix", Dim = c(M,M),
        ## 'i :' maybe there's a faster way (w/o matrix indexing), but elegant?
        i = as.vector(matrix(0L:(M-1L), nrow=k)[, rep(seq_len(N), each=k)]),
        p = k * 0L:M,
        x = as.double(unlist(lmat, recursive=FALSE, use.names=FALSE)))
  }
  
  ESSphi1 <- matrix(0,n,K)
  ESSphi2 <- array(0,c(K,K,n))
  
  ESSpsi1 <- matrix(0,np,K)
  ESSpsi2 <- array(0,c(K,K,np))
  
  ESSphipsi2 <- array(0,c(K,K,np))
  
  
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
      N <- m[j]/nperiods[j]           #Assumes equal number of subjects/period...maybe sum over cluster-periods?????
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
      obs = c(t(residj))
      tm1 <- nperiods[j]*(N-1)*log(det(SigmaE))+(nperiods[j]-1)*log(det(SigmaE+N*SigmaPsi))+log(det(SigmaE+N*(SigmaPsi+nperiods[j]*SigmaPhi)))
      InvSS2 <- solve(SigmaE+N*SigmaPsi)-InvS2E
      InvSS22 <- solve(SigmaE+N*(SigmaPsi+nperiods[j]*SigmaPhi))-solve(SigmaE+N*SigmaPsi)
      Invj <- kronecker(diag(1,nrow=nperiods[j]),kronecker(diag(1,nrow=N),InvS2E) + 
        kronecker(matrix(1,N,N),InvSS2)/N) +
        kronecker(matrix(1,nperiods[j],nperiods[j]),kronecker(matrix(1,N,N),InvSS22)/(nperiods[j]*N))
      tm2 <- c(t(obs) %*% Invj %*% obs)
      temp <- temp-(tm1+tm2)/2
    }
    
    temp
  }
  thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaPsi[!lower.tri(SigmaPsi)]),c(SigmaE[!lower.tri(SigmaE)]))
  LLold <- loglik(thetah)
  
  
  niter=1
  while((niter <= maxiter) & (abs(delta) > epsilon)){
    
    # Expectation step
    count <- 1
    count2 <- 1
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      Zj <- Z[ID == j,] 
      Zj <- Zj[,colSums(Zj)>0] 
      Zjj <- matrix(0,2*nrow(Zj),2*ncol(Zj))
      for (k in 1:ncol(Zj)){
        for (r in 1:nrow(Zj)){
          Zjj[r,2*k-1] <- Zj[r,k]
          Zjj[nrow(Zj)+r,2*k] <- Zj[r,k]
        }
      }
      mlist <- c(replicate(1,InvS2Phi,simplify=FALSE),replicate(nperiods[j], InvS2Psi, simplify=FALSE))
      InvS2Zeta <- bdiag_m(mlist)
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
      #ZZt <- crossprod(Zj)
      #Vj <- solve(InvS2Zeta + m[j]*kronecker(InvS2E, ZZt))
      Vj <- solve(InvS2Zeta + m[j]*t(Zjj)%*%kronecker(diag(1,nrow=nrow(Xj)),InvS2E)%*%Zjj)
      r1 <- t(Zj)%*%residj[,1]
      r2 <- t(Zj)%*%residj[,2]
      Muj <- Vj %*% rbind(InvS2E[1,1]*r1 + InvS2E[1,2]*r2,
                              InvS2E[2,1]*r1 + InvS2E[2,2]*r2)
      Nujj <- Vj + tcrossprod(Muj)
      
      ESSphi1[j,] <- Muj[1:K,]
      ESSphi2[,,j] <- as.matrix(Nujj[1:K,1:K])
      
      Sij <- Muj[-c(1:K),]
      VSij <- Nujj[,-(1:K)][-(1:K),]
      
      for(k in 1:nperiods[j]){
        ESSpsi1[count,] <- Sij[1:K]
        ESSpsi2[,,count] <- as.matrix(VSij[1:K,1:K])
        count <- count + 1
        Sij <- Sij[-(1:K)]
        VSij <- VSij[,-(1:K)][-(1:K),]
      }
      
      Vbs <- Nujj
      for(k in 1:t){
        ESSphipsi2[,,count2] <- as.matrix(Vbs[1:K,(K+1):(2*K)]) + as.matrix(Vbs[(K+1):(2*K),1:K]) 
        Vbs <- Vbs[,-((K+1):(2*K))][-((K+1):(2*K)),]
        count2 <- count2 + 1
      }
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
      rowSums(sweep(ESSphipsi2,3,mp,FUN="*"),dims=2)
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
    if(verbose) cat(paste('iter=',niter),'\t',paste('loglik=',LLnew),'\t',paste('SigmaE=',c(SigmaE[!lower.tri(SigmaE)])),'\n')
    
    #print(niter)
    #print(zeta)
    #print(SigmaPhi)
    #print(SigmaE)
    #print(LLnew)
  }
  param <- list(theta=list(zeta=zeta,SigmaE=SigmaE,SigmaPhi=SigmaPhi,SigmaPsi=SigmaPsi),loglik=LLnew,eps=epsilon,iter=niter)
  return(param) 
}







