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
  
  vc1<-as.data.frame(VarCorr(fm1))
  vc2<-as.data.frame(VarCorr(fm2))
  s2e1 <- ifelse(vc1[vc1$grp=="Residual",4]==0,0.01,vc1[vc1$grp=="Residual",4])
  s2e2 <- ifelse(vc2[vc2$grp=="Residual",4]==0,0.01,vc2[vc2$grp=="Residual",4])
  SigmaE <- diag(c(s2e1, s2e2))
  InvS2E <- solve(SigmaE)
  
  fe1 = as.numeric(fixef(fm1))
  fe2 = as.numeric(fixef(fm2))
  delta1<-fe1[length(fe1)]
  delta2<-fe2[length(fe2)]
  delta<-(delta1/sqrt(s2e1) + delta2/sqrt(s2e2))/2                                  # Better way to initialize?
  zeta<-c(fe1[-length(fe1)],fe2[-length(fe2)],delta)
  beta1<-c(fe1[-length(fe1)],delta)
  beta2<-c(fe2[-length(fe2)],delta)
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
  
  s2phi1 <- ifelse(vc1[vc1$grp==paste(cluster),4]==0,0.01,vc1[vc1$grp==paste(cluster),4])
  s2phi2 <- ifelse(vc2[vc2$grp==paste(cluster),4]==0,0.01,vc2[vc2$grp==paste(cluster),4])
  SigmaPhi <- diag(c(s2phi1, s2phi2))
  InvS2Phi <- solve(SigmaPhi)
  
  s2psi1 <- ifelse(vc1[vc1$grp==paste(cluster.period),4]==0,0.01,vc1[vc1$grp==paste(cluster.period),4])
  s2psi2 <- ifelse(vc2[vc2$grp==paste(cluster.period),4]==0,0.01,vc2[vc2$grp==paste(cluster.period),4])
  SigmaPsi <- diag(c(s2psi1, s2psi2))
  InvS2Psi <- solve(SigmaPsi)
  
  
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
    beta1 = c(theta[1:(nvar-1)],theta[(2*nvar-1)])
    beta2 = theta[nvar:(2*nvar-1)]
    sphi11 = theta[(2*nvar)]
    sphi12 = theta[(2*nvar+1)]
    sphi22 = theta[(2*nvar+2)]
    se11 = theta[(2*nvar+3)]
    se12 = theta[(2*nvar+4)]
    se22 = theta[(2*nvar+5)]
    spsi11 = theta[(2*nvar+6)]
    spsi12 = theta[(2*nvar+7)]
    spsi22 = theta[(2*nvar+8)]
    SigmaPhi = matrix(c(sphi11,sphi12,sphi12,sphi22),2,2)
    SigmaPsi = matrix(c(spsi11,spsi12,spsi12,spsi22),2,2)
    SigmaE = matrix(c(se11,se12,se12,se22),2,2)
    InvS2Phi <- solve(SigmaPhi)
    InvS2Psi <- solve(SigmaPsi)
    InvS2E <- solve(SigmaE)
    
    temp <- 0
    for(j in 1:n){
      Yj <- Y[ID == j,,drop=FALSE]
      Xj <- X[ID == j,,drop=FALSE]
      residj <- Yj - cbind(sqrt(se11)*Xj%*%beta1, sqrt(se22)*Xj%*%beta2)
      obs = c(t(residj))
      psizes <- mp[cID.period==j]
      if(m[j]/nperiods[j] == psizes[1]){     #If equal number of subjects/period
        N <- psizes[1]          
        tm1 <- nperiods[j]*(N-1)*log(det(SigmaE))+(nperiods[j]-1)*log(det(SigmaE+N*SigmaPsi))+log(det(SigmaE+N*(SigmaPsi+nperiods[j]*SigmaPhi)))
        InvSS2 <- solve(SigmaE+N*SigmaPsi)-InvS2E
        InvSS22 <- solve(SigmaE+N*(SigmaPsi+nperiods[j]*SigmaPhi))-solve(SigmaE+N*SigmaPsi)
        Invj <- kronecker(diag(1,nrow=nperiods[j]),kronecker(diag(1,nrow=N),InvS2E) + 
                            kronecker(matrix(1,N,N),InvSS2)/N) +
          kronecker(matrix(1,nperiods[j],nperiods[j]),kronecker(matrix(1,N,N),InvSS22)/(nperiods[j]*N))
        tm2 <- c(t(obs) %*% Invj %*% obs)
        temp <- temp-(tm1+tm2[1])/2
      } else {                                #If unequal number of subjects/period
        bdiag_s <- matrix(0,nrow=2*sum(psizes),ncol=2*sum(psizes))
        last.row<-0
        for (k in 1:nperiods[j]){
          bdiag_s[(last.row+1):(last.row+2*psizes[k]),(last.row+1):(last.row+2*psizes[k])] <- kronecker(matrix(1,psizes[k],psizes[k]), SigmaPsi)
          last.row<-last.row+2*psizes[k]
        }
        Omega <- bdiag_s + kronecker(diag(1,m[j]),SigmaE) + kronecker(matrix(1,m[j],m[j]),SigmaPhi)
        Invj <- solve(Omega)
        tm1 <- log(det(Omega))
        tm2 <- t(obs) %*% Invj %*% obs
        temp <- temp-(tm1+tm2[1])/2
      }
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
      #Zjj <- kronecker(Zj,diag(1,nrow=K))
      Zjj <- kronecker(Zj,diag(c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2])),nrow=K))
      mlist <- c(replicate(1,InvS2Phi,simplify=FALSE),replicate(nperiods[j], InvS2Psi, simplify=FALSE))
      InvS2Zeta <- bdiag_m(mlist)
      residj <- Yj - cbind(sqrt(SigmaE[1,1])*Xj%*%beta1, sqrt(SigmaE[2,2])*Xj%*%beta2)
      residjj <- kronecker(residj[,1],matrix(c(1,0),nrow=K)) + kronecker(residj[,2],matrix(c(0,1),nrow=K))
      Vj <- solve(InvS2Zeta + t(Zjj)%*%kronecker(diag(1,nrow=nrow(Xj)),InvS2E)%*%Zjj)
      Muj <- Vj %*% t(Zjj)%*%kronecker(diag(1,nrow=nrow(Xj)),InvS2E)%*%residjj
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
    #X.fe <- cbind(kronecker(diag(c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2])),2,2),X[,-nvar]),c(sqrt(SigmaE[1,1])*X[,nvar],sqrt(SigmaE[2,2])*X[,nvar]))        #KP: new design matrix in order to get one delta!!!!
    X.fe <- cbind(kronecker(X[,-nvar],diag(c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2])),nrow=K)),kronecker(X[,nvar],c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2]))))
    mlist <- replicate(nrow(X.fe)/K,InvS2E,simplify=FALSE)
    InvS2Zeta <- bdiag_m(mlist)
    Vzeta <- solve(t(X.fe)%*%InvS2Zeta%*%X.fe)
    rzeta1<-Y[,1]-sqrt(SigmaE[1,1])*ESSphi1[ID,1]-sqrt(SigmaE[1,1])*ESSpsi1[ID.period,1]
    rzeta2<-Y[,2]-sqrt(SigmaE[2,2])*ESSphi1[ID,2]-sqrt(SigmaE[2,2])*ESSpsi1[ID.period,2]
    rzeta <- kronecker(rzeta1,matrix(c(1,0),nrow=K)) + kronecker(rzeta2,matrix(c(0,1),nrow=K))
    zeta <- Vzeta %*% (t(X.fe)%*%InvS2Zeta%*%rzeta)
    zeta <- as.vector(zeta)
    beta1 = zeta[c(TRUE,FALSE)]
    beta2 = c(zeta[c(FALSE,TRUE)],zeta[length(zeta)])
    zeta <- c(beta1[-length(beta1)],beta2)
    
    # Maximization step - epsilon
    re <- Y - cbind(sqrt(SigmaE[1,1])*X%*%beta1, sqrt(SigmaE[2,2])*X%*%beta2)
    rss <- crossprod(re) + rowSums(sweep(ESSphi2,3,m,FUN="*"),dims=2)*matrix(c(SigmaE[1,1],sqrt(SigmaE[1,1])*sqrt(SigmaE[2,2]),sqrt(SigmaE[1,1])*sqrt(SigmaE[2,2]),SigmaE[2,2]),nrow=2) + 
      rowSums(sweep(ESSpsi2,3,mp,FUN="*"),dims=2)*matrix(c(SigmaE[1,1],sqrt(SigmaE[1,1])*sqrt(SigmaE[2,2]),sqrt(SigmaE[1,1])*sqrt(SigmaE[2,2]),SigmaE[2,2]),nrow=2) -
      crossprod(sweep(ESSphi1,2,c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2])),FUN="*"),rowsum(re,ID)) - crossprod(rowsum(re,ID),sweep(ESSphi1,2,c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2])),FUN="*")) -
      crossprod(sweep(ESSpsi1,2,c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2])),FUN="*"),rowsum(re,ID.period)) - crossprod(rowsum(re,ID.period),sweep(ESSpsi1,2,c(sqrt(SigmaE[1,1]),sqrt(SigmaE[2,2])),FUN="*")) +
      rowSums(sweep(ESSphipsi2,3,mp,FUN="*"),dims=2)*matrix(c(SigmaE[1,1],sqrt(SigmaE[1,1])*sqrt(SigmaE[2,2]),sqrt(SigmaE[1,1])*sqrt(SigmaE[2,2]),SigmaE[2,2]),nrow=2)
    SigmaE <- rss/sum(m)
    InvS2E <- solve(SigmaE)
    
    # whether the algorithm converges
    thetah = c(zeta,c(SigmaPhi[!lower.tri(SigmaPhi)]),c(SigmaE[!lower.tri(SigmaE)]),c(SigmaPsi[!lower.tri(SigmaPsi)]))
    LLnew <- loglik(thetah)
    delta <- abs(LLnew - LLold)
    LLold <- LLnew
    converge = (abs(delta)<=epsilon)
    niter <- niter + 1
    if(verbose) cat(paste('iter=',niter),'\t',paste('param.error=',epsilon),'\t',paste('loglik=',LLnew),'\n');
    #if(verbose) cat(paste('iter=',niter),'\t',paste('loglik=',LLnew),'\t',paste('SigmaE=',c(SigmaE[!lower.tri(SigmaE)])),'\n')
    
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
  
  param <- list(theta=list(zeta=zeta,SigmaE=SigmaE,SigmaPhi=SigmaPhi,SigmaPsi=SigmaPsi),loglik=LLnew,eps=epsilon,iter=niter,
                Vtheta=Vtheta,SEcheck=SEcheck) #SEtheta=SEtheta,
  return(param) 
}







