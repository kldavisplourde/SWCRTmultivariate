# Application Study

source("/Users/kdavis07/Documents/GitHub/CoPrimarySWCRT/powerSampleCal_IU_HoopGir.R")

# Parameter Inputs
sd1<-sqrt(611.13)
sd2<-sqrt(695.73)
rho2<-matrix(c(1,0.58,0.58,1),2)
rho02<-matrix(c(0.006,0,0,0.029),2);rho01<-matrix(c(0.00002,0,0,0.0068),2);deltas<-c(0.3*sd1,0.35*sd2);t<-5;N=(t-1)*4;m<-12

# Power
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
# I=16 and N=12 gives 86.3% power

##### Sensitivity Analysis for Application Study #####
### CAC=0.2
rho01<-matrix(c(0.00122,0,0,0.0059),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### CAC=0.5
rho01<-matrix(c(0.003,0,0,0.015),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{12} = -0.002
rho01<-matrix(c(0.00122,0,0,0.0059),2)
rho02<-matrix(c(0.006,-0.002,-0.002,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{12} = 0.002
rho02<-matrix(c(0.006,0.002,0.002,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{2} = 0.023
rho02<-matrix(c(0.006,0,0,0.023),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{2} = 0.035
rho02<-matrix(c(0.006,0,0,0.035),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.005
rho02<-matrix(c(0.005,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.007
rho02<-matrix(c(0.007,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.46
rho02<-matrix(c(0.006,0,0,0.029),2)
rho2<-matrix(c(1,0.46,0.46,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.7
rho2<-matrix(c(1,0.7,0.7,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)

#Additional scenarios (up to 40%)
### CAC=0
rho01<-matrix(c(0,0,0,0),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### CAC=0.8
rho01<-matrix(c(0.0049,0,0,0.0235),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{12} = -0.004
rho01<-matrix(c(0.00122,0,0,0.0059),2)
rho02<-matrix(c(0.006,-0.004,-0.004,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{12} = 0.004
rho02<-matrix(c(0.006,0.004,0.004,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{2} = 0.017
rho02<-matrix(c(0.006,0,0,0.017),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{2} = 0.041
rho02<-matrix(c(0.006,0,0,0.041),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.004
rho02<-matrix(c(0.004,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.008
rho02<-matrix(c(0.008,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.35
rho02<-matrix(c(0.006,0,0,0.029),2)
rho2<-matrix(c(1,0.35,0.35,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.81
rho2<-matrix(c(1,0.81,0.81,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)

#Additional scenarios x2 (up to 60%)
### rho_0^{2} = 0.012
rho02<-matrix(c(0.006,0,0,0.012),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{2} = 0.046
rho02<-matrix(c(0.006,0,0,0.046),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.002
rho02<-matrix(c(0.002,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_0^{1} = 0.010
rho02<-matrix(c(0.010,0,0,0.029),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.23
rho02<-matrix(c(0.006,0,0,0.029),2)
rho2<-matrix(c(1,0.23,0.23,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)
### rho_2^{12} = 0.93
rho2<-matrix(c(1,0.93,0.93,1),2)
calPower_IU(deltas,margins=c(0,0),vars=c(sd1^2,sd2^2),rho01,rho02,rho2,N,t,m,K=2,alpha=0.05)