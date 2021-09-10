library(nlme)
library(foreach)
library(doMC)
library(doRNG)
#library(doSNOW)

#load.lib<-c("rmutil","dplyr", "foreach", "reshape2","doSNOW",  "doParallel","mvtnorm","nlme","numDeriv")
#install.lib<-load.lib[!load.lib %in% installed.packages()]
#for(lib in install.lib) install.packages(lib,dependencies=TRUE ,repo="http://cran.rstudio.com/")
#sapply(load.lib, require, character=TRUE)

setwd("/Users/kdavis07/Documents/GitHub/CoPrimarySWCRT")
source("gendata_copri_varCluster_HH.R")
source("EM_uncorrected_HH.R")

args<-commandArgs(trailingOnly = TRUE)
k<-as.integer(args[1])
if (is.na(k)) k <- 1
paste("Scenario:",k)

ncores<-as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK",1))
if (is.na(ncores)) ncores<-1
registerDoMC(cores=ncores)

# define scenarios
scenarios <- read.table("/Users/kdavis07/Dropbox/SW-CRT Methods Development/2_CoPrimary/RCode/Simulations/Preliminary/prelim_params.txt", header=TRUE, sep="")
scenarios <- subset(scenarios, scenario == k)

scenario <- k
t <- scenarios$t
N <- scenarios$N
cs <- scenarios$m
eff<-c(scenarios$delta1,scenarios$delta2)
rho01<-matrix(c(scenarios$rho01.11,scenarios$rho01.12,scenarios$rho01.12,scenarios$rho01.22),2)
rho2<-matrix(c(1,scenarios$rho2.12,scenarios$rho2.12,1),2)
vars<-c(scenarios$var1,scenarios$var2)
bs <- 0
beta <- cumsum(c(0.1,0.1*0.5,0.1*(0.5^2),0.1*(0.5^3),0.1*(0.5^4),0.1*(0.5^5)))[1:t-1]
nsim<-2


for(a in 1:1)
  #foreach(a=1,.combine=rbind, .verbose = T,.export=c('datagen_cont', 'EM.estim'), .packages=c( "mvtnorm","numDeriv","nlme"))%do%
{
  set.seed(k*1838+t*10+N+cs)
  fail_count <- 0
  max_fail <- 200
  
  ZETA <- SIGMAE <- SIGMAPHI <- SEtheta <- results <- NULL
  # Loop Index
  i<-0
  # While Loop (discarding false data)
  while(i<nsim){
    i<-i+1
    itemp<-i
    set.seed(i+173)
  
    dt<-datagen_cont(n=N, m=cs, K=2, cv=0, rho01=rho01, rho2=rho2, vars=vars, eff=eff, time.eff=c(beta,beta))
    data<-dt$short

    if(t==3){
      lme1<-lme(out1~time.1+time.2+arm,random=~1|cluster,data=data,control=lmeControl(returnObject=TRUE))
      lme2<-lme(out2~time.1+time.2+arm,random=~1|cluster,data=data,control=lmeControl(returnObject=TRUE))
    }
  
    if(t==4){
      lme1<-lme(out1~time.1+time.2+time.3+arm,random=~1|cluster,data=data,control=lmeControl(returnObject=TRUE))
      lme2<-lme(out2~time.1+time.2+time.3+arm,random=~1|cluster,data=data,control=lmeControl(returnObject=TRUE))
    }
  
    if(t==5){
      lme1<-lme(out1~time.1+time.2+time.3+time.4+arm,random=~1|cluster,data=data,control=lmeControl(returnObject=TRUE))
      lme2<-lme(out2~time.1+time.2+time.3+time.4+arm,random=~1|cluster,data=data,control=lmeControl(returnObject=TRUE))
    }

    param<-try(EM.estim(data,lme1,lme2, maxiter=500, epsilon=1e-4, verbose=FALSE))
    if(class(param)=="try-error"){i<- i-1;fail_count <-fail_count+1}
    if(i<itemp){next}
    if(fail_count > max_fail){break}
    ZETA <- rbind(ZETA,param$theta$zeta)
    SIGMAE[[i]] <- param$theta$SigmaE
    SIGMAPHI[[i]] <- param$theta$SigmaPhi
    
    thetah<- c(param$theta$zeta,c(param$theta$SigmaPhi[!lower.tri(param$theta$SigmaPhi)]),param$theta$SigmaE[!lower.tri(param$theta$SigmaE)])

    TermsX1 <- lme1$terms
    TermsX2 <- lme2$terms
    mfX1 <- model.frame(TermsX1, data = data)[,-1]
    mfX2 <- model.frame(TermsX2, data = data)[,-1]
    Y <- as.matrix(cbind(model.frame(TermsX1, data = data)[,1],model.frame(TermsX2, data = data)[,1]))
    X <- as.matrix(cbind(1, mfX1)) # design matrix
    
    ID <- lme1$groups[[1]]
    n <- length(unique(ID))
    m <- as.numeric(table(lme1$groups[[1]]))
    nvar<-length(as.numeric(lme1$coefficients$fixed))
    
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
    
    Vtheta = try(solve(-hessian(loglik,thetah)))
    if(class(Vtheta)=="try-error"){i<- i-1; fail_count <- fail_count+1}
    if(i<itemp){next}
    if(fail_count > max_fail){break}
    SEtheta = rbind(SEtheta, sqrt(diag(Vtheta)))
  
  #c(betas,SigmaPhi,SigmaE,SEtheta,covDelta12,covDelta21)
  results<-rbind(results,c(thetah,SEtheta))
  }
  
  simData <- results
}
  #stopCluster(makeCluster(ncores))

if(t==3){
  colnames(simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se")
}

if(t==4){
  colnames(simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se")
}

if(t==5){
  colnames(simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Period5.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Period5.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Period5.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Period5.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se")
}

if(sum(eff)==0) analysis<-"error"
if(sum(eff) != 0) analysis<-"power"

simData <- as.data.frame(simData)
write.table(simData, file=paste("results/UncorrectedResults_",analysis,"_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)
#write.table(simData, file=paste("results/CorrectedResults_",analysis,"_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)