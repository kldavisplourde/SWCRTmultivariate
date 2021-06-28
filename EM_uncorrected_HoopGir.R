#setwd("C:/Users/lifan/Dropbox/CRT_Coprimary/Data/")
#data=read.csv("short2.csv",header=TRUE)

#library(nlme)
library(mvtnorm)
library(numDeriv)

# function to perform EM estimation with K=2 outcomes
EM.estim <- function(data, fm1,fm2, maxiter=500,epsilon=1e-4
                     , verbose=FALSE){
  # fit mixed model to initialize parameters
  #fm1 <- lme(formula1, random = ~ 1|cluster, data=data,control=lmeControl(returnObject=TRUE))
  #fm2 <- lme(formula2, random = ~ 1|cluster, data=data,control=lmeControl(returnObject=TRUE))
  K <- 2
  zeta <- as.numeric(c(fm1$coefficients$fixed, fm2$coefficients$fixed))
  beta1 = as.numeric(fm1$coefficients$fixed)
  beta2 = as.numeric(fm2$coefficients$fixed)
  if (length(as.numeric(fm1$coefficients$fixed)) != length(as.numeric(fm2$coefficients$fixed)))
    stop("\nnumber of covariates do not match between endpoints.")
  nvar<-length(as.numeric(fm1$coefficients$fixed))
  TermsX1 <- fm1$terms
  TermsX2 <- fm2$terms
  mfX1 <- model.frame(TermsX1, data = data)[,-1]
  mfX2 <- model.frame(TermsX2, data = data)[,-1]
  if (identical(mfX1,mfX2) == FALSE)
    stop("\ncovariates do not match between endpoints.")
  ##
  
  # vector of cluster sizes
  m <- as.numeric(table(fm1$groups[[1]]))
  
  s2phi1 <- VarCorr(fm1)[1,1]
  s2phi2 <- VarCorr(fm2)[1,1]
  SigmaPhi <- diag(c(s2phi1, s2phi2)) # KP: Initialized assuming independence
  InvS2Phi <- solve(SigmaPhi)
  
  s2e1 <- VarCorr(fm1)[2,1]
  s2e2 <- VarCorr(fm2)[2,1]
  SigmaE <- diag(c(s2e1, s2e2))
  InvS2E <- solve(SigmaE)
  
  Y <- as.matrix(cbind(model.frame(TermsX1, data = dt)[,1],model.frame(TermsX2, data = dt)[,1]))
  ID <- fm1$groups[[1]]
  n <- length(unique(ID))
  X <- as.matrix(cbind(1, mfX1)) # design matrix
  
  ESSphi1 <- matrix(0,n,K)
  ESSphi2 <- array(0,c(K,K,n))
  
  
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
    SigmaPhi = matrix(c(sphi11,sphi12,sphi12,sphi22),2,2)
    SigmaE = matrix(c(se11,se12,se12,se22),2,2)
    InvS2Phi <- solve(SigmaPhi)
    InvS2E <- solve(SigmaE)
    
    temp <- 0
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
      obs = c(t(residj))
      tm1 <- (m[j]-1)*log(det(SigmaE))+log(det(SigmaE+m[j]*SigmaPhi))
      InvSS2 <- solve(SigmaE+m[j]*SigmaPhi)-InvS2E
      Invj <- kronecker(diag(nrow=m[j]),InvS2E) + 
        kronecker(matrix(1,m[j],m[j]),InvSS2)/m[j]
      tm2 <- c(t(obs) %*% Invj %*% obs)
      temp <- temp-(tm1+tm2)/2
    }
    temp
  }
  thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]))
  LLold <- loglik(thetah)
  
  
  niter=1
  while((niter <= maxiter) & (abs(delta) > epsilon)){
    
    # Expectation step
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(Xj%*%beta1, Xj%*%beta2)
      Vj <- solve(InvS2Phi + m[j]*InvS2E)
      Muj <- as.numeric(Vj %*% InvS2E %*% colSums(residj))
      Nujj <- Vj + tcrossprod(Muj)
      ESSphi1[j,] <- Muj
      ESSphi2[,,j] <- Nujj
    }
    
    # Maximization step - phi
    SigmaPhi <- apply(ESSphi2,1:2, sum)/n
    InvS2Phi <- solve(SigmaPhi)
    
    # Maximization step - zeta
    # Simplify the expression analytically, and obtain simple expression!
    XXt <- crossprod(X)
    Vzeta <-solve(kronecker(InvS2E, XXt))
    rzeta1 <- t(X)%*%(Y[,1]-ESSphi1[ID,1])
    rzeta2 <- t(X)%*%(Y[,2]-ESSphi1[ID,2])
    zeta <- Vzeta %*% rbind(InvS2E[1,1]*rzeta1 + InvS2E[1,2]*rzeta2,
                            InvS2E[2,1]*rzeta1 + InvS2E[2,2]*rzeta2)
    zeta <- c(zeta)
    beta1 = zeta[1:nvar]
    beta2 = zeta[(nvar+1):(2*nvar)]
    
    # Maximization step - epsilon
    re <- Y - cbind(X%*%beta1, X%*%beta2)
    rss <- crossprod(re) + rowSums(sweep(ESSphi2,3,m,FUN="*"),dims=2) -
      crossprod(ESSphi1,rowsum(re,ID)) - crossprod(rowsum(re,ID),ESSphi1)
    SigmaE <- rss/sum(m)
    # SigmaE <- diag(diag(SigmaE))
    InvS2E <- solve(SigmaE)
    
    # whether the algorithm converges
    thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]))
    LLnew <- loglik(thetah)
    delta <- abs(LLnew - LLold)
    LLold <- LLnew
    converge = (abs(delta)<=epsilon)
    niter <- niter + 1
    if(verbose) cat(paste('iter=',niter),'\t',paste('param.error=',epsilon),'\t',paste('loglik=',LLnew),'\n');  
    
    #print(niter)
    #print(zeta)
    #print(SigmaPhi)
    #print(SigmaE)
    #print(LLnew)
  }
  param <- list(theta=list(zeta=zeta,SigmaE=SigmaE,SigmaPhi=SigmaPhi),loglik=LLnew,eps=epsilon,iter=niter)
  return(param) 
}







