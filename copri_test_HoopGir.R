library(lme4)
#library(nlme)
library(doMC)
library(doRNG)
library(lmeInfo)

setwd("/Users/kdavis07/Documents/GitHub/CoPrimarySWCRT")
source("gendata_copri_varCluster_HoopGir.R")
source("EM_uncorrected_HoopGir.R")

args<-commandArgs(trailingOnly = TRUE)
k<-as.integer(args[1])
if (is.na(k)) k <- 1
paste("Scenario:",k)

ncores<-as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK",1))
if (is.na(ncores)) ncores<-1
registerDoMC(cores=ncores)

# define scenarios
scenarios <- read.table("/Users/kdavis07/Dropbox/SW-CRT Methods Development/2_CoPrimary/RCode/Simulations/HoopGir/Sim_Params.txt", header=TRUE, sep="")
scenarios <- subset(scenarios, scenario == k)

scenario <- k
t <- scenarios$t
N <- scenarios$N
cs <- scenarios$m
eff<-c(scenarios$delta1,scenarios$delta2)
rho01<-matrix(c(scenarios$rho01.11,scenarios$rho01.12,scenarios$rho01.12,scenarios$rho01.22),2)
rho02<-matrix(c(scenarios$rho02.11,scenarios$rho02.12,scenarios$rho02.12,scenarios$rho02.22),2)
rho2<-matrix(c(1,scenarios$rho2.12,scenarios$rho2.12,1),2)
vars<-c(4,4) #c(scenarios$var1,scenarios$var2)
bs <- 0
beta <- cumsum(c(0.1,0.1*0.5,0.1*(0.5^2),0.1*(0.5^3),0.1*(0.5^4),0.1*(0.5^5)))[1:t-1]
nsim<-2

set.seed(3826+k)
fail_count <- 0
max_fail <- 200

simData <- naive.simData <- NULL
# Loop Index
i<-0
# While Loop (discarding false data)
while(i<nsim){
  i<-i+1
  itemp<-i
  
  dt<-datagen_cont(n=N, m=cs, K=2, cv=0, rho01=rho01, rho02=rho02, rho2=rho2, vars=vars, eff=eff, time.eff=c(beta,beta))
  data<-dt$short

  if(t==3){
    lme1<-lmer(out1~time.1+time.2+arm +(1|cluster) +(1|cluster.period),data=data)
    lme2<-lmer(out2~time.1+time.2+arm +(1|cluster) +(1|cluster.period),data=data)
  }
  
  if(t==4){
    lme1<-lmer(out1~time.1+time.2+time.3+arm +(1|cluster) +(1|cluster.period),data=data)
    lme2<-lmer(out2~time.1+time.2+time.3+arm +(1|cluster) +(1|cluster.period),data=data)
  }
  
  if(t==5){
    lme1<-lmer(out1~time.1+time.2+time.3+time.4+arm +(1|cluster) +(1|cluster.period),data=data)
    lme2<-lmer(out2~time.1+time.2+time.3+time.4+arm +(1|cluster) +(1|cluster.period),data=data)
  }
  
  param<-try(EM.estim(data,lme1,lme2,cluster="cluster",cluster.period="cluster.period",maxiter=500, epsilon=1e-4, verbose=FALSE))
  if(class(param)=="try-error"){
    i<- i-1;fail_count <-fail_count+1
  }else{
    if(param$SEcheck=="ERROR"|anyNA(param$theta$zeta)==TRUE|anyNA(param$theta$SigmaE)==TRUE|anyNA(param$theta$SigmaPhi)==TRUE|anyNA(param$theta$SigmaPsi)==TRUE){i<- i-1;fail_count <-fail_count+1}
  }
  if(fail_count > max_fail){break}
  if(i<itemp){next}
  
  SEtheta.i<-sqrt(diag(param$Vtheta))
  results.i<- c(param$theta$zeta,param$theta$SigmaPhi[!lower.tri(param$theta$SigmaPhi)],param$theta$SigmaPsi[!lower.tri(param$theta$SigmaPsi)],param$theta$SigmaE[!lower.tri(param$theta$SigmaE)],SEtheta.i)
  
  # From individual models (not taking into account between-outcome within-subject correlation)
  naive.zeta <- as.numeric(c(fixef(lme1), fixef(lme2)))
  
  vc1<-as.data.frame(VarCorr(lme1))
  vc2<-as.data.frame(VarCorr(lme2))
  naive.SigmaPhi <- c(vc1[vc1$grp=="cluster",4],vc2[vc2$grp=="cluster",4])
  naive.SigmaPsi <- c(vc1[vc1$grp=="cluster.period",4],vc2[vc2$grp=="cluster.period",4])
  naive.SigmaE <- c(vc1[vc1$grp=="Residual",4],vc2[vc2$grp=="Residual",4])
  
  naive.SE <- as.numeric(c(sqrt(diag(vcov(lme1))),sqrt(diag(vcov(lme2))))) # add SE for random effects?

  naive.i<-c(naive.zeta,naive.SigmaPhi,naive.SigmaPsi,naive.SigmaE,naive.SE)
  
  # combining results
  simData<-rbind(simData,results.i)
  naive.simData<-rbind(naive.simData,naive.i)
}

if(t==3){
  colnames(simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22",
                       "SigmaPsi11","SigmaPsi12","SigmaPsi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se",
                       "SigmaPsi11.se","SigmaPsi12.se","SigmaPsi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se")
  
  colnames(naive.simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Treatment.est1",
                             "Intercept.est2","Period2.est2","Period3.est2","Treatment.est2",
                             "SigmaPhi11","SigmaPhi22","SigmaPsi11","SigmaPsi22","SigmaE11","SigmaE22",
                             "Intercept.se1","Period2.se1","Period3.se1","Treatment.se1",
                             "Intercept.se2","Period2.se2","Period3.se2","Treatment.se2")
}

if(t==4){
  colnames(simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22",
                       "SigmaPsi11","SigmaPsi12","SigmaPsi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se",
                       "SigmaPsi11.se","SigmaPsi12.se","SigmaPsi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se")
  
  colnames(naive.simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Treatment.est1",
                             "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Treatment.est2",
                             "SigmaPhi11","SigmaPhi22","SigmaPsi11","SigmaPsi22","SigmaE11","SigmaE22",
                             "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Treatment.se1",
                             "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Treatment.se2")
}

if(t==5){
  colnames(simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Period5.est1","Treatment.est1",
                       "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Period5.est2","Treatment.est2",
                       "SigmaPhi11","SigmaPhi12","SigmaPhi22",
                       "SigmaPsi11","SigmaPsi12","SigmaPsi22","SigmaE11","SigmaE12","SigmaE22",
                       "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Period5.se1","Treatment.se1",
                       "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Period5.se2","Treatment.se2",
                       "SigmaPhi11.se","SigmaPhi12.se","SigmaPhi22.se",
                       "SigmaPsi11.se","SigmaPsi12.se","SigmaPsi22.se","SigmaE11.se","SigmaE12.se","SigmaE22.se")
  
  colnames(naive.simData)<-c("Intercept.est1","Period2.est1","Period3.est1","Period4.est1","Period5.est1","Treatment.est1",
                             "Intercept.est2","Period2.est2","Period3.est2","Period4.est2","Period5.est2","Treatment.est2",
                             "SigmaPhi11","SigmaPhi22","SigmaPsi11","SigmaPsi22","SigmaE11","SigmaE22",
                             "Intercept.se1","Period2.se1","Period3.se1","Period4.se1","Period5.se1","Treatment.se1",
                             "Intercept.se2","Period2.se2","Period3.se2","Period4.se2","Period5.se2","Treatment.se2")
}

if(0 %in% eff) analysis<-"error" else analysis<-"power"

simData <- as.data.frame(simData)
naive.simData <- as.data.frame(naive.simData)

write.table(simData, file=paste("results/UncorrectedResults_",analysis,"_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)
write.table(naive.simData, file=paste("results/NaiveResults_",analysis,"_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)
#write.table(simData, file=paste("results/CorrectedResults_",analysis,"_scenario",scenario,".txt",sep=""), sep="\t", row.names=F)