
# This code is used for simulating the extended MDS model described in Liu et al., 2020, 
# Computational Brain & Behavior. 
# Compiling MDS.c is required before running this code. 
# This code can recreate figure 5,6,7,8 for simulation 1 in the paper
# A previously simulated data at line 81 can be used for recreating the figures without simulating the data
# Questions can be sent to Qingfang Liu, psychliuqf@gmail.com
# Feb, 2020

rm(list = ls())

## MDS function complied from MDS.c 
MDS=function(pL,n_region,n_time,n_trial,n_rep,thresh1,thresh2,sigma1,sigma2,
             c0,c1,c2,a1,a2,a3,d1,d2,tau,L,ds_factor,A12,A34,A5,A6,more=FALSE){
  
  dyn.load("MDS.so")
  out=.C('MDS_two_choice',pL=as.double(pL),n_region=as.integer(n_region),n_time=as.integer(n_time),
         n_trial=as.integer(n_trial),n_rep=as.integer(n_rep),thresh1=as.double(thresh1),thresh2=as.double(thresh2),
         sigma1=as.double(sigma1),sigma2=as.double(sigma2),c0=as.double(c0),c1=as.double(c1),c2=as.double(c2),
         a1=as.double(a1),a2=as.double(a2),a3=as.double(a3),d1=as.double(d1),d2=as.double(d2),
         tau=as.double(tau),L=as.integer(L),ds_factor=as.integer(ds_factor),
         A12=as.double(A12),A34=as.double(A34),A5=as.double(A5),A6=as.double(A6),
         Resp=numeric(n_trial*n_rep),t1=numeric(n_trial*n_rep),t0=numeric(n_trial*n_rep),
         y_BOLD=numeric(n_rep*n_time*n_trial*n_region/ds_factor)
         )
  dat=list("Resp"=array(out$Resp,dim = c(n_rep,n_trial)),
       "RT"=array(out$t1+tau,dim = c(n_rep,n_trial)),
       "BOLD"=array(out$y_BOLD,c(n_time*n_trial/ds_factor,n_region,n_rep))
       )
  
  if(more==TRUE){ # print t0 and t1 for visulization
    dat[["t0"]]=array(out$t0,dim = c(n_rep,n_trial))
    dat[["t1"]]=array(out$t1,dim = c(n_rep,n_trial))
  }
  dat
}   

##################### simulation setup ############################3

conds = seq(.1,.9,.1) # nine PL conditions
n_trial_per_cond = 30 # number of trials in each condition

true=NULL # true parameter values for simulation in the paper
set.seed(3)
true$pL=sample(rep(conds,n_trial_per_cond)) # sequence of interleaved conditions
true$n_region=6 # number of regions
true$n_time=2000 # number of time points
true$n_trial=length(true$pL) # number of trials in each experiment
true$n_rep=2 # number of replication of one experiment

true$thresh1=250 # theta_1 from the paper
true$thresh2=1500 # theta_2 from the paper
true$sigma1=16 
true$sigma2=5 
true$c0=.7  
true$c1=.5  
true$c2=.9 
true$a1=.8 
true$a2=-.2 
true$a3=-.8 
true$d1=.9 
true$d2=.9 
true$tau=100 
true$A12=.0005 
true$A34=.00006 
true$A5=.0015 
true$A6=.0002 
true$L=32000 # length of HRF in msec unit (i.e. 32s)
true$ds_factor=1000 # downsample factor to approx real fMRI temporal resolution

############################ simulating data #####################

dat=MDS(true$pL,true$n_region,true$n_time,true$n_trial,true$n_rep,true$thresh1,true$thresh2,true$sigma1,true$sigma2,
        true$c0,true$c1,true$c2,true$a1,true$a2,true$a3,true$d1,true$d2,true$tau,true$L,true$ds_factor,
        true$A12,true$A34,true$A5,true$A6,more = TRUE)

save.image(file = "MDS_simu.RData")

################## Analyzing simulated data ##################################

load(file = "MDS_simu.RData")

# aggregate data according to conditions

Resp_by_cond = sapply(conds,function(x) dat$Resp[,(true$pL)==x])
Resp_by_cond[Resp_by_cond==9]=NA
RT_by_cond = sapply(conds,function(x) dat$RT[,(true$pL)==x])
RT_by_cond[is.na(Resp_by_cond)]=NA

# Fig.5 - choice response time distributions

layout(matrix(c(1,1,1,2,3,6,9,12,4,7,10,12,5,8,11,12),nrow = 4,byrow = F),widths = c(.8,3,3,3),heights = c(3,3,3,1))
mar1=c(.1,.1,.1,.1)
mar2=c(2,2.5,2,1)
cex.big=3
cex.tiny=2

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,6),xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")
text(.8,3,"Density",xlim=c(0,1),ylim=c(0,6),cex=cex.big,srt=90)

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")

par(mar=mar2)
for(i in c(1:length(conds))){
  RT=RT_by_cond[,i]
  choice=Resp_by_cond[,i]
  RT<-RT[!is.na(RT)]
  choice=choice[!is.na(choice)]
  
  RT[choice==1]=RT[choice==1]*(-1)
  hist(RT,breaks = seq(-2100,2100,by = 100),col = "light blue",xaxs="r",yaxs="r",freq = F,
       xlab = "",cex=2,cex.axis=1.8,
       cex.lab=1.5,main = '',xlim = c(-2100,2100),ylim = c(0,.007))
  if(i==1)
  {text(900,.005,expression(paste(p[L]," = .1")),cex = 3)
    text(-1000,.003,"L",cex = 3)
    text(1000,.003,"R",cex = 3)
  }
  if(i==2)text(900,.005,expression(paste(p[L]," = .2")),cex = 3)
  if(i==3)text(900,.005,expression(paste(p[L]," = .3")),cex = 3)
  if(i==4)text(900,.005,expression(paste(p[L]," = .4")),cex = 3)
  if(i==5)text(900,.005,expression(paste(p[L]," = .5")),cex = 3)
  if(i==6)text(900,.005,expression(paste(p[L]," = .6")),cex = 3)
  if(i==7)text(900,.005,expression(paste(p[L]," = .7")),cex = 3)
  if(i==8)text(900,.005,expression(paste(p[L]," = .8")),cex = 3)
  if(i==9)text(900,.005,expression(paste(p[L]," = .9")),cex = 3)
}

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,6),xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")
text(.5,3,"Response Times (msec)",xlim=c(0,1),ylim=c(0,6),cex=cex.big)

# Fig. 6 - mean response times and accuracy by conditions

mean_RT_by_cond = apply(RT_by_cond,2,mean,na.rm=T)
names(mean_RT_by_cond) = seq(.1,.9,.1)
acc = apply(Resp_by_cond,2,mean,na.rm=T)-1
acc[6:9]=1-acc[6:9]
std <- function(x) sd(x,na.rm = T)/sqrt(length(x[!is.na(x)]))
std_RT_by_cond = apply(RT_by_cond,2,std)

par(mfrow=c(1,2),mar=c(5,5,2,2))
plot(seq(.1,.9,.1),acc,pch=16,cex=2,ylab = "Accuracy",xlab = expression(paste(p[L])),cex.lab=1.5)
ymax = max(mean_RT_by_cond)+3*std_RT_by_cond[which.max(mean_RT_by_cond)]

barCenters <- barplot(height = mean_RT_by_cond,names.arg = seq(.1,.9,.1),
                      beside = true, las = 1,ylim = c(0, ymax),cex.names = 0.95,
                      xlab = expression(paste(p[L])),ylab = "Mean Response Time (msec)",
                      border = "black", axes = TRUE,col = "light blue",cex.lab=1.5)
segments(barCenters, mean_RT_by_cond - std_RT_by_cond * 2, barCenters,
         mean_RT_by_cond + std_RT_by_cond * 2, lwd = 1.5)
arrows(barCenters, mean_RT_by_cond - std_RT_by_cond * 2, barCenters,
       mean_RT_by_cond + std_RT_by_cond * 2, lwd = 1.5, angle = 90,
       code = 3, length = 0.05)


# Fig. 7 -  BOLD activity across regions 

# colors for plotting
col.n = "orange"
col.r = c("#CCFFCC","#FFCC66","#66CC00","#FF9900","#333366","#FF3300")
library(plotrix)

set.seed(2)
rep=sample(true$n_rep,1)

layout(matrix(c(1,2,3,1,4,5,1,6,7,8,9,9),nrow = 4,byrow = TRUE),widths = c(.5,3,3),heights = c(3,3,3,1))
mar1=c(.8,1.5,.5,1.2)
mar2=c(.1,.1,.3,.1)
cex.big=3
cex.tiny=2
cex.little=1.8

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,6),xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")
text(.5,3,"Simulated BOLD Signal from MDS",xlim=c(0,1),ylim=c(0,6),cex=cex.big,srt=90)

par(mar=mar1)
plot(dat$BOLD[,1,rep],type = "l",xlab = "",main ="" ,ylab = "BOLD response",cex.main=cex.tiny,cex.axis=cex.little,xaxt="n")
points(dat$BOLD[,1,rep],col=2,pch=16)
draw.circle(510,1.55,20,col = col.r[1])
text(510,1.55,expression(R[1]),cex = 2.5)

par(mar=mar1)
plot(dat$BOLD[,2,rep],type = "l",xlab = "",main = "",ylab = "BOLD response",cex.main=cex.tiny,cex.axis=cex.little,xaxt="n")
points(dat$BOLD[,2,rep],col=2,pch=16)
draw.circle(510,1.65,20,col = col.r[2])
text(510,1.65,expression(R[2]),cex = 2.5)

par(mar=mar1)
plot(dat$BOLD[,3,rep],type = "l",xlab = "",main = "",ylab = "BOLD response",cex.main=cex.tiny,cex.axis=cex.little,xaxt="n")
points(dat$BOLD[,3,rep],col=2,pch=16)
draw.circle(510,1.62,20,col = col.r[3])
text(510,1.62,expression(R[3]),cex = 2.5)

par(mar=mar1)
plot(dat$BOLD[,4,rep],type = "l",xlab = "",main = "",ylab = "BOLD response",cex.main=cex.tiny,cex.axis=cex.little,xaxt="n")
points(dat$BOLD[,4,rep],col=2,pch=16)
draw.circle(510,1.64,20,col = col.r[4])
text(510,1.64,expression(R[4]),cex = 2.5)

par(mar=mar1)
plot(dat$BOLD[,5,rep],type = "l",xlab = "",main = "",ylab = "BOLD response",cex.main=cex.tiny,cex.axis=cex.little)
points(dat$BOLD[,5,rep],col=2,pch=16)
draw.circle(510,1.35,20,col = col.r[5])
text(510,1.35,expression(R[5]),cex = 2.5)

par(mar=mar1)
plot(dat$BOLD[,6,rep],type = "l",xlab = "",main = "",ylab = "BOLD response",cex.main=cex.tiny,cex.axis=cex.little)
points(dat$BOLD[,6,rep],col=2,pch=16)
draw.circle(510,.85,20,col = col.r[6])
text(510,.85,expression(R[6]),cex = 2.5)

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,6),xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")
text(.5,3,"Time (sec)",xlim=c(0,1),ylim=c(0,6),cex=cex.big)


# Fig 8 - GLM fit for simulated BOLD signal

library(plotrix)

rep_cond = c(sapply(true$pL,function(x) rep(x,2)))

cond_BOLD=NULL
mean_BOLD = matrix(NA,9,6)
for(i in 1:9){
  cond_BOLD[[i]] = dat$BOLD[rep_cond==.1*i,,]
  mean_BOLD[i,] = apply(cond_BOLD[[i]],2,mean)
}

hrf = function(ts,A=1){ # HRF convolution function
  a1 = 6; a2 = 16; b1 = 1; b2 = 1; c = 1/6
  ys = A * ( ((ts^(a1-1) * b1^a1 * exp(-b1 * ts) )/gamma(a1) )
             - c * (((ts^(a2 - 1) * b2^a2 * exp(-b2 * ts) )/gamma(a2))) )
  ys
}

L=32
y_hrf= hrf(seq(1, L, 1))
y_hrf = c(y_hrf,rep(0,length(rep_cond)-1)) # zero pad
tmp_cond = c(rep_cond,rep(0,L-1))*10

conv_sig = Re(fft(fft(tmp_cond)*fft(y_hrf),inverse = T))[1:length(rep_cond)]
conv_sig = conv_sig/(length(rep_cond)+L-1)


layout(matrix(c(1,1,2,3,6,9,4,7,9,5,8,9),3,4),widths = c(.5,2,2,2),heights = c(2,2,.5))
ylab="Simulated BOLD Signal from MDS Model (Y)"

mar1=c(2,2,1,1)
cex.big=3
cex.tiny=2.5
cex.little=1.5

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,6),xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")
text(.5,3,ylab,xlim=c(0,1),ylim=c(0,6),cex=cex.big,srt=90)

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,6),xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")

par(mar=mar1)
plot(conv_sig,apply(dat$BOLD[,1,],1,mean),col=col.n,cex.axis=cex.little)
draw.circle(6.5,.85,.7,col = col.r[1])
text(6.5,.85,expression(R[1]),cex = 2.5)
abline(lm(apply(dat$BOLD[,1,],1,mean) ~ conv_sig),col=2)
summary(lm(apply(dat$BOLD[,1,],1,mean) ~ conv_sig))

par(mar=mar1)
plot(conv_sig,apply(dat$BOLD[,2,],1,mean),col=col.n,cex.axis=cex.little)
draw.circle(6.5,.85,.7,col = col.r[2])
text(6.5,.85,expression(R[2]),cex = 2.5)
abline(lm(apply(dat$BOLD[,2,],1,mean) ~ conv_sig),col=2)
summary(lm(apply(dat$BOLD[,2,],1,mean) ~ conv_sig))

par(mar=mar1)
plot(conv_sig,apply(dat$BOLD[,3,],1,mean),col=col.n,cex.axis=cex.little)
draw.circle(6.5,.79,.7,col = col.r[3])
text(6.5,.79,expression(R[3]),cex = 2.5)
abline(lm(apply(dat$BOLD[,3,],1,mean) ~ conv_sig),col=2)
summary(lm(apply(dat$BOLD[,3,],1,mean) ~ conv_sig))

par(mar=mar1)
plot(conv_sig,apply(dat$BOLD[,4,],1,mean),col=col.n,cex.axis=cex.little)
draw.circle(6.5,.78,.7,col = col.r[4])
text(6.5,.78,expression(R[4]),cex = 2.5)
abline(lm(apply(dat$BOLD[,4,],1,mean) ~ conv_sig),col=2)
summary(lm(apply(dat$BOLD[,4,],1,mean) ~ conv_sig))

par(mar=mar1)
plot(conv_sig,apply(dat$BOLD[,5,],1,mean),col=col.n,cex.axis=cex.little)
draw.circle(6.5,-.05,.7,col = col.r[5])
text(6.5,-.05,expression(R[5]),cex = 2.5)
abline(lm(apply(dat$BOLD[,5,],1,mean) ~ conv_sig),col=2)
summary(lm(apply(dat$BOLD[,5,],1,mean) ~ conv_sig))

par(mar=mar1)
plot(conv_sig,apply(dat$BOLD[,6,],1,mean),col=col.n,cex.axis=cex.little)
draw.circle(6.5,.22,.7,col = col.r[6])
text(6.5,.22,expression(R[6]),cex = 2.5)
abline(lm(apply(dat$BOLD[,6,],1,mean) ~ conv_sig),col=2)
summary(lm(apply(dat$BOLD[,6,],1,mean) ~ conv_sig))

par(mar=mar1)
plot(NA,xlim=c(0,1),ylim=c(0,6),xlab="",ylab="",main="",xaxt="n",yaxt="n",bty="n")
text(.5,3,expression("(convolved) "*p[L]*" (X)"),xlim=c(0,1),ylim=c(0,6),cex=cex.big)

