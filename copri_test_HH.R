library(nlme)
library(doMC)
library(doRNG)
library(lmeInfo)

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
scenarios <- read.table("/Users/kdavis07/Dropbox/SW-CRT Methods Development/2_CoPrimary/RCode/Simulations/Sim_Params.txt", header=TRUE, sep="")
scenarios <- subset(scenarios, scenario == k)

scenario <- k
t <- scenarios$t
N <- scenarios$N
cs <- scenarios$m
eff<-c(scenarios$delta1,scenarios$delta2)
rho01<-matrix(c(scenarios$rho01.11,scenarios$rho01.12,scenarios$rho01.12,scenarios$rho01.22),2)
rho2<-matrix(c(1,scenarios$rho2.12,scenarios$rho2.12,1),2)
vars<-c(1,1) #c(scenarios$var1,scenarios$var2)
bs <- 0
beta <- cumsum(c(0.1,0.1*0.5,0.1*(0.5^2),0.1*(0.5^3),0.1*(0.5^4),0.1*(0.5^5)))[1:t-1]
nsim<-2

set.seed(8374+k)
fail_count <- 0
max_fail <- 200
  
simData <- naive.simData <- NULL
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
  if(param$SEcheck=="ERROR"|anyNA(param$theta$zeta)==TRUE|anyNA(param$theta$SigmaE)==TRUE|anyNA(param$theta$SigmaPhi)==TRUE){i<- i-1;fail_count <-fail_count+1}
  if(fail_count > max_fail){break}
  if(i<itemp){next}
  
  results.i<- c(param$theta$zeta,c(param$theta$SigmaPhi[!lower.tri(param$theta$SigmaPhi)]),param$theta$SigmaE[!lower.tri(param$theta$SigmaE)],param$SEtheta)
  
  # From individual models (not taking into account between-outcome within-subject correlation)
  naive.zeta <- as.numeric(c(lme1$coefficients$fixed, lme2$coefficients$fixed))
  naive.SE <- as.numeric(c(sqrt(diag(lme1$varFix)),sqrt(diag(lme2$varFix)),sqrt(diag(varcomp_vcov(lme1))),sqrt(diag(varcomp_vcov(lme2)))))
  naive.SigmaE<-as.numeric(c(VarCorr(lme1)[2,1],VarCorr(lme2)[2,1]))
  naive.SigmaPhi<-as.numeric(c(VarCorr(lme1)[1,1],VarCorr(lme2)[1,1]))
  
  naive.i<-c(naive.zeta,naive.SigmaPhi,naive.SigmaE,naive.SE)
  
  # combining results
  simData<-rbind(simData,results.i)
  naive.simData<-rbind(naive.simData,naive.i)
}
#stopCluster(makeCluster(ncores))

if(t==3){
  colnames(simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se")
  
  colnames(naive.simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi22","SigmaE11","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaE11.se","SigmaPhi22.se","SigmaE22.se")
}

if(t==4){
  colnames(simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se")
  
  colnames(naive.simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi22","SigmaE11","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaE11.se","SigmaPhi22.se","SigmaE22.se")
}

if(t==5){
  colnames(simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Period5.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Period5.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Period5.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Period5.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se")
  
  colnames(naive.simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Period5.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Period5.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi22","SigmaE11","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Period5.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Period5.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaE11.se","SigmaPhi22.se","SigmaE22.se")
}

if(sum(eff)==0) analysis<-"error"
if(sum(eff) != 0) analysis<-"power"

simData <- as.data.frame(simData)
naive.simData <- as.data.frame(naive.simData)

write.table(simData, file=paste("results/UncorrectedResults_",analysis,"_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)
write.table(naive.simData, file=paste("results/NaiveResults_",analysis,"_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)
#write.table(simData, file=paste("results/CorrectedResults_",analysis,"_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)