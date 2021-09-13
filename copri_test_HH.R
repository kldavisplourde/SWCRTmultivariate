library(nlme)
#library(foreach)
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

set.seed(k*697+t*4)
fail_count <- 0
max_fail <- 200
  
ZETA <- SIGMAE <- SIGMAPHI <- SEtheta <- simData <- NULL
# Loop Index
i<-0
# While Loop (discarding false data)
while(i<nsim){
  i<-i+1
  itemp<-i
  #set.seed(i+173)
  
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
  if(anyNA(param$theta$zeta)==TRUE|anyNA(param$theta$SigmaE)==TRUE|anyNA(param$theta$SigmaPhi)==TRUE|anyNA(param$SEtheta)==TRUE){i<- i-1;fail_count <-fail_count+1}
  if(fail_count > max_fail){break}
  if(i<itemp){next}
  ZETA <- rbind(ZETA,param$theta$zeta)
  SIGMAE[[i]] <- param$theta$SigmaE
  SIGMAPHI[[i]] <- param$theta$SigmaPhi
  SEtheta <- rbind(SEtheta,param$SEtheta)
    
  results.i<- c(param$theta$zeta,c(param$theta$SigmaPhi[!lower.tri(param$theta$SigmaPhi)]),param$theta$SigmaE[!lower.tri(param$theta$SigmaE)],SEtheta)
  
  simData<-rbind(simData,results.i)
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