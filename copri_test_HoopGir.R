library(nlme)
library(foreach)
#set.seed(3628) #1000 iterations
set.seed(5792) #10000 iterations

simData <- foreach(i=1:10000, .combine=rbind) %do% {
  dt<-datagen_cont(n=15, m=5, K=2, cv=0, sigmac=matrix(c(1.3,0.15,0.15,1.1),2), sigmae= matrix(c(2,0.1,0.1,1.5),2), eff=c(1,0.5), time.eff=c(1.1,1.3,1.4,0.9,0.7,0.5))$short

  lme1<-lme(out1~time.1+time.2+time.3+arm,random=~1|cluster,data=dt,control=lmeControl(returnObject=TRUE))
  lme2<-lme(out2~time.1+time.2+time.3+arm,random=~1|cluster,data=dt,control=lmeControl(returnObject=TRUE))
  #formula1=formula(lme1)
  #formula2=formula(lme2)

  fitEM<-EM.estim(dt,lme1,lme2)

  betas<-fitEM$theta$zeta
  SigmaE<-c(fitEM$theta$SigmaE[!lower.tri(fitEM$theta$SigmaE)])
  SigmaPhi<-c(fitEM$theta$SigmaPhi[!lower.tri(fitEM$theta$SigmaPhi)])
  iter<-fitEM$iter
  
  c(betas,SigmaPhi,SigmaE,iter)
}

colnames(simData)<-c("Intercept.est1","Time1.est1","Time2.est1","Time3.est1","Treatment.est1",
                     "Intercept.est2","Time1.est2","Time2.est2","Time3.est2","Treatment.est2",
                     "SigmaPhi11","SigmaPhi12","SigmaPhi22","SigmaE11","SigmaE12","SigmaE22","EM.iter")

simData <- as.data.frame(simData)
write.table(simData, file="/Users/kdavis07/Desktop/test.EM.results10k2.txt", sep="\t", row.names=F)