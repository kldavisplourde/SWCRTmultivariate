# Original code from Fan Li, modified by Kendra Plourde

###############################################
### function for generating the SW-CRT data ###
###############################################
#n=15; m=5; K=2; cv=0; sigmac=matrix(c(1,0.01,0.01,1),K); sigmacp=matrix(c(1,0.01,0.01,1),K); sigmae= matrix(c(1,0.01,0.01,1),K); eff=c(1,0.5); time.eff=c(1.1,1.3,1.4,0.9,0.7,0.5)
datagen_cont <- function(n, m, K, cv, sigmac, sigmacp, sigmae, eff, time.eff ){
  
  # n -number of clusters
  # m - cluster size
  # K - number of multiple co-primary outcomes
  # cv - cluster size coefficient of variation
  # sigmac - correlation matrix of random cluster effects
  # sigmacp - correlation matrix of random cluster-period effects
  # sigmae - correlation matrix of random errors
  # eff - a vector of treatment effect sizes 
  # time.eff - vector of time effects (order is by primary outcome and period)
  library(MASS)
  
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
  #cluster <- rep(1:(2*Narm),  mj);
  cluster <- rep(1:n,  mj*nperiods);
  
  
  #generate patient id
  patid <- seq(1:n)
  temp <- cbind(patid, cluster)
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
  mjp <- rep(mj, each=nperiods)                             # KP: Assuming cluster size does not change over time.
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
  
 