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
  
  s2phi1 <- VarCorr(fm1)[1,1]
  s2phi2 <- VarCorr(fm2)[1,1]
  SigmaPhi <- diag(c(s2phi1, s2phi2)) # KP: Initialized assuming independence
  InvS2Phi <- solve(SigmaPhi)
  
  s2e1 <- as.numeric(VarCorr(fm1)[2,1])
  s2e2 <- as.numeric(VarCorr(fm2)[2,1])
  SigmaE <- diag(c(s2e1, s2e2))
  InvS2E <- solve(SigmaE)
  
  zeta <- as.numeric(c(fm1$coefficients$fixed, fm2$coefficients$fixed))
  ## KP: We now have a time effect to take into account which depends on the number of periods.
  #beta1 = zeta[1:2] 
  #beta2 = zeta[3:4]
  fe1 = as.numeric(fm1$coefficients$fixed)
  fe2 = as.numeric(fm2$coefficients$fixed)
  delta1<-fe1[length(fe1)]
  delta2<-fe2[length(fe2)]
  delta<-(delta1/sqrt(s2e1) + delta2/sqrt(s2e2))/2                                  # Better way to initialize?
  zeta<-c(fe1[-length(fe1)],fe2[-length(fe2)],delta)
  beta1<-c(fe1[-length(fe1)],delta)
  beta2<-c(fe2[-length(fe2)],delta)
  
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
  #m <- as.numeric(table(data$cluster))
  m <- as.numeric(table(fm1$groups[[1]]))
  
  #Y <- as.matrix(data[,c("out1","out2")])
  Y <- as.matrix(cbind(model.frame(TermsX1, data = data)[,1],model.frame(TermsX2, data = data)[,1]))
  #ID <- as.numeric(data$cluster)
  ID <- fm1$groups[[1]]
  n <- length(unique(ID))
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
  
  
  #maxiter=500
  #epsilon=1e-4
  delta = 2*epsilon
  max_modi = 20
  
  converge = 0
  
  # log likelihood
  
  loglik = function(theta){
    beta1 = c(theta[1:(nvar-1)],theta[(2*nvar-1)])
    beta2 = theta[nvar:(2*nvar-1)]
    sphi11 = theta[(2*nvar)]
    sphi12 = theta[(2*nvar+1)]
    sphi22 = theta[(2*nvar+2)]
    se11 = theta[(2*nvar+3)]
    se12 = theta[(2*nvar+4)]
    se22 = theta[(2*nvar+5)]
    SigmaPhi = matrix(c(sphi11,sphi12,sphi12,sphi22),2,2)
    SigmaE = matrix(c(se11,se12,se12,se22),2,2)
    InvS2Phi <- solve(SigmaPhi)
    InvS2E <- solve(SigmaE)
    
    temp <- 0
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(sqrt(se11)*Xj%*%beta1, sqrt(se22)*Xj%*%beta2)
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
      residj <- Yj - cbind(sqrt(SigmaE[1,1])*Xj%*%beta1, sqrt(SigmaE[2,2])*Xj%*%beta2)
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
    X.fe <- cbind(kronecker(X[,-nvar],diag(c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2])),nrow=K)),kronecker(X[,nvar],c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2]))))
    mlist <- replicate(nrow(X.fe)/K,InvS2E,simplify=FALSE)
    InvS2Zeta <- bdiag_m(mlist)
    Vzeta <- solve(t(X.fe)%*%InvS2Zeta%*%X.fe)
    rzeta1<-Y[,1]-sqrt(SigmaE[1,1])*ESSphi1[ID,1]
    rzeta2<-Y[,2]-sqrt(SigmaE[2,2])*ESSphi1[ID,2]
    rzeta <- kronecker(rzeta1,matrix(c(1,0),nrow=K)) + kronecker(rzeta2,matrix(c(0,1),nrow=K))
    zeta <- Vzeta %*% (t(X.fe)%*%InvS2Zeta%*%rzeta)
    zeta <- as.vector(zeta)
    beta1 = zeta[c(TRUE,FALSE)]
    beta2 = c(zeta[c(FALSE,TRUE)],zeta[length(zeta)])
    zeta <- c(beta1[-length(beta1)],beta2)
    
    # Maximization step - epsilon
    re <- Y - cbind(sqrt(SigmaE[1,1])*X%*%beta1, sqrt(SigmaE[2,2])*X%*%beta2)
    rss <- crossprod(re) + rowSums(sweep(ESSphi2,3,m,FUN="*"),dims=2)*matrix(c(SigmaE[1,1],sqrt(SigmaE[1,1])*sqrt(SigmaE[2,2]),sqrt(SigmaE[1,1])*sqrt(SigmaE[2,2]),SigmaE[2,2]),nrow=2) -
      crossprod(sweep(ESSphi1,2,c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2])),FUN="*"),rowsum(re,ID)) - crossprod(rowsum(re,ID),sweep(ESSphi1,2,c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2])),FUN="*"))
    SigmaE <- rss/sum(m)
    # SigmaE <- diag(diag(SigmaE))
    InvS2E <- solve(SigmaE)
    
    # whether the algorithm converges
    # thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),diag(SigmaE))
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
  
  Vtheta <- matrix(NA,length(thetah),length(thetah))
  try(Vtheta <- solve(-hessian(loglik,thetah))) #, silent = TRUE
  SEcheck <- "GOOD"
  if(anyNA(Vtheta)==TRUE|min(diag(Vtheta))<0) SEcheck <- "ERROR"
  #Vtheta = try(solve(-hessian(loglik,thetah)))
  #SEtheta = sqrt(diag(Vtheta))
  #if(class(Vtheta)=="try-error"){i<- i-1; fail_count <- fail_count+1}
  #if(i<itemp){next}
  #if(fail_count > max_fail){break}
  #SEtheta = rbind(SEtheta, sqrt(diag(Vtheta)))
  
  param <- list(theta=list(zeta=zeta,SigmaE=SigmaE,SigmaPhi=SigmaPhi),loglik=LLnew,eps=epsilon,iter=niter,
                Vtheta=Vtheta,SEcheck=SEcheck) #SEtheta=SEtheta,
  return(param) 
}







