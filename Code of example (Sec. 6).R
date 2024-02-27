rm(list=ls())

install.packages("foreach")
install.packages("parallel")
install.packages("doParallel")

library(foreach)
library(parallel)
library(doParallel)

cl=makeCluster(detectCores(logical = F))
registerDoParallel(cl)

Total.Data=read.csv("New_York_City.csv")
New_York.Data=Total.Data[which(Total.Data$County=="New York City"),]
start=which(New_York.Data$Test.Date=="03/01/2022")
end=which(New_York.Data$Test.Date=="07/15/2022")

#-------------Test.Data----------------
Test.Data=New_York.Data[end:start,c(1,3,5,7)]
F_Test.D=Test.Data
for(i in 1:(dim(Test.Data)[1])){
  F_Test.D[i,]=Test.Data[dim(Test.Data)[1]+1-i,]
}

t1=which(F_Test.D$Test.Date=="03/01/2022")
t2=which(F_Test.D$Test.Date=="04/01/2022")
t3=which(F_Test.D$Test.Date=="05/01/2022")
t4=which(F_Test.D$Test.Date=="06/01/2022")
t5=which(F_Test.D$Test.Date=="07/01/2022")
xais.l=c(t1,t2,t3,t4,t5)
Phase.t=which(F_Test.D$Test.Date=="03/20/2022")


##------------Plot------------------------------
cT=c(1:dim(F_Test.D)[1])
plot(cT,F_Test.D$Test...Positive,type="l",xlab="Time",ylab="Positive rate",
xaxt='n',lwd=0.7,cex.main=0.8)
axis(1,at=xais.l,labels=c("03/01","04/01","05/01","06/01","07/01"))
abline(v=Phase.t,lty=3,col="black",lwd=2)
ss=lm(F_Test.D$Test...Positive[Phase.t+1:dim(F_Test.D)[1]]~cT[Phase.t+1:dim(F_Test.D)[1]])
summary(ss)
curve(0.0093467+0.0005489*x,from=Phase.t+1,add=T,lty=4,col="red")

par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(2,1,0))
plot(cT,F_Test.D$Total.Number.of.Tests.Performed,type="l",xlab="Time",ylab="Total number of tests (1,000)",
     xaxt='n',main="a",lwd=0.7,cex.main=0.8)
axis(1,at=xais.l,labels=c("03/01","04/01","05/01","06/01","07/01"))
c1= 69.47 ;c2=150.19 ;c3= -41.92
curve(c1/(1+exp(-(x-c2)/c3)),add=T,lty=4,col="red")
lm(F_Test.D$Total.Number.of.Tests.Performed~cT)
#curve( 72.0731-0.1997*x,add=T,lty=4,col="red")
plot(cT,F_Test.D$New.Positives,type="l",xlab="Time",ylab="New positives (1,000)",
     xaxt='n',main="b",lwd=0.7,cex.main=0.8)
axis(1,at=xais.l,labels=c("03/01","04/01","05/01","06/01","07/01"))
abline(v=Phase.t,lty=3,col="black",lwd=2)


#--------model selection------------
# linear regression or logistic model
Sample.size=F_Test.D$Total.Number.of.Tests.Performed
year=c(1:length(Sample.size))
nls.data=data.frame(cbind(Sample.size-1,year))  #model:nt=logistic+1
pop.mod1=nls(Sample.size ~ SSlogis(year,phi1,phi2,phi3),data=nls.data)
#59212/length(F_Test.D$New.Positives)
library(ggplot2)
p <- ggplot(nls.data,aes(year, Sample.size))
p+geom_point(size=3)+geom_line(aes(year,fitted(pop.mod1)),col='red')
pop.mod2=lm(Sample.size~year)
mean(resid(pop.mod2)^2)
num=F_Test.D$New.Positives[1:Phase.t]
mean(num);var(num)



#-------Function_Test----------------
lambda0=mean(F_Test.D$Test...Positive[1:Phase.t])
lambda1=max(F_Test.D$Test...Positive[1:Phase.t])
c_theta.D=c(0.001,0.0025)
gamma=0.002
ARL0=1/0.0027;
h.low=0;h.up=10;
N=20000

#------------------Function-------------------------
CUSUM.S.hf=function(lambda1,lambda0,ARL0,N,h.low,h.high){
  
  c1= 69.47 ;c2=150.19 ;c3= -41.92
  nt.f=function(t){
    nt=c1/(1+exp(-(t-c2)/c3))
    return(nt)
  }
  
  CSM.h=function(lambda1,lambda0){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      t=t+1
      n[t]=nt.f(t)
      x[t]=rpois(1,lambda0*n[t])
      if(t==1){W[t]=max(0,x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))
      } else{W[t]=max(0,W[t-1]+x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))}
      plot.stat=W[t]
    }
    t
  }
  
  for(i in 1:100){
    
    h=(h.high+h.low)/2
    RL=foreach(i=1:N,.combine="c")%dopar% CSM.h(lambda1,lambda0)
    
    
    if(abs(mean(RL)-ARL0)<2 | h.high-h.low<0.00005){
      return(h)
      break
    }else if(mean(RL)-ARL0>2){
      h.high=h
    }else{h.low=h}
    print(h)
  }
  return(h)
}
CUSUM.D.hf=function(theta.D,lambda0,ARL0,N,h.low,h.high){
  
  c1= 69.47 ;c2=150.19 ;c3= -41.92
  nt.f=function(t){
    nt=c1/(1+exp(-(t-c2)/c3))
    return(nt)
  }
  
  CSD.h=function(theta.D,lambda0){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      t=t+1
      n[t]=nt.f(t)
      x[t]=rpois(1,lambda0*n[t])
      lambda_star=lambda0+t*theta.D*sqrt(lambda0)
      
      if(t==1){W[t]=max(0,x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))
      } else{W[t]=max(0,W[t-1]+x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))}
      plot.stat=W[t]
    }
    t
  }
  
  for(i in 1:100){
    
    h=(h.high+h.low)/2
    RL=foreach(i=1:N,.combine="c")%dopar% CSD.h(theta.D,lambda0)
    
    
    if(abs(mean(RL)-ARL0)<2 | h.high-h.low<0.00005){
      return(h)
      break
    }else if(mean(RL)-ARL0>2){
      h.high=h
    }else{h.low=h}
    print(h)
  }
  return(h)
}
ACUSUM_VN1.h.f=function(gamma,lambda0,ARL0,N,h.low,h.high){
  c1= 69.47 ;c2=150.19 ;c3= -41.92
  nt.f=function(t){
    nt=c1/(1+exp(-(t-c2)/c3))
    return(nt)
  }
  ACSM.1=function(gamma,lambda0){
    
    W=n=x=y1=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(t)
      x[t]=rpois(1,lambda0*n[t])
      
      if(t==1){
        y1[t]=max(lambda0,(1-gamma)*lambda0+gamma*x[t]/n[t])
        W[t]=max(0,x[t]*log(y1[t]/lambda0)-n[t]*(y1[t]-lambda0))
        
        
      } else{
        y1[t]=max(lambda0,(1-gamma)* y1[t-1]+gamma*x[t]/n[t])
        W[t]=max(0,W[t-1]+x[t]*log(y1[t]/lambda0)-n[t]*(y1[t]-lambda0))
        
      }
      plot.stat=W[t]
    }
    
    t
  }
  
  
  for(i in 1:100){
    
    h=(h.high+h.low)/2
    RL=foreach(i=1:N,.combine="c")%dopar%ACSM.1(gamma,lambda0)
    
    
    if(abs(mean(RL)-ARL0)<2 | h.high-h.low<0.00005){
      return(h)
      break
    }else if(mean(RL)-ARL0>2){
      h.high=h
    }else{h.low=h}
    print(c(h.high,h.low))
  }
  return(h)
}
ACUSUM_VN2.h.f=function(gamma,lambda0,ARL0,N,h.low,h.high){
  
  c1= 69.47 ;c2=150.19 ;c3= -41.92
  nt.f=function(t){
    nt=c1/(1+exp(-(t-c2)/c3))
    return(nt)
  }
  ASM.2=function(gamma,lambda0){
    
    W=n=x=y1=y2=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(t)
      x[t]=rpois(1,lambda0*n[t])
      
      if(t==1){
        y1[t]=(1-gamma)*lambda0+gamma*x[t]/n[t]
        y2[t]=max(lambda0,y1[t])
        W[t]=max(0,x[t]*log(y2[t]/lambda0)-n[t]*(y2[t]-lambda0))
      } else{
        y1[t]=(1-gamma)*y1[t-1]+gamma*x[t]/n[t]
        y2[t]=max(lambda0,y1[t])
        W[t]=max(0,W[t-1]+x[t]*log(y2[t]/lambda0)-n[t]*(y2[t]-lambda0))
      }
      plot.stat=W[t]
    }
    t
  }
  
  for(i in 1:100){
    
    h=(h.high+h.low)/2
    RL=foreach(i=1:N,.combine="c")%dopar% ASM.2(gamma,lambda0)
    
    
    if(abs(mean(RL)-ARL0)<2 | h.high-h.low<0.00005){
      return(h)
      break
    }else if(mean(RL)-ARL0>2){
      h.high=h
    }else{h.low=h}
    print(c(h.high,h.low))
  }
  return(h)
}
ACUSUM_VN3.h.f=function(gamma,lambda0,ARL0,N,h.low,h.high){
  
  c1= 69.47 ;c2=150.19 ;c3= -41.92
  nt.f=function(t){
    nt=c1/(1+exp(-(t-c2)/c3))
    return(nt)
  }
  
  ASM.3=function(gamma,lambda0){
    
    W=n=x=y3=y4=y_xt=y_nt=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(t)
      x[t]=rpois(1,lambda0*n[t])
      
      if(t==1){
        y_xt[t]=gamma*x[t]+(1-gamma)*lambda0*n[t]
        y_nt[t]=gamma*n[t]+(1-gamma)*n[t]
        y3[t]=y_xt[t]/y_nt[t]
        y4[t]=max(lambda0,y3[t])
        W[t]=max(0,x[t]*log(y4[t]/lambda0)-n[t]*(y4[t]-lambda0))
      } else{
        y_xt[t]=gamma*x[t]+(1-gamma)*y_xt[t-1]
        y_nt[t]=gamma*n[t]+(1-gamma)*y_nt[t-1]
        y3[t]=y_xt[t]/y_nt[t]
        y4[t]=max(lambda0,y3[t])
        W[t]=max(0,W[t-1]+x[t]*log(y4[t]/lambda0)-n[t]*(y4[t]-lambda0))
      }
      plot.stat=W[t]
    }
    t
  }
  
  for(i in 1:100){
    
    h=(h.high+h.low)/2
    RL=foreach(i=1:N,.combine="c")%dopar%ASM.3(gamma,lambda0) 
    
    
    if(abs(mean(RL)-ARL0)<2 | h.high-h.low<0.00005){
      return(h)
      break
    }else if(mean(RL)-ARL0>2){
      h.high=h
    }else{h.low=h}
    print(c(h.high,h.low))
  }
  return(h)
}

H_s=CUSUM.S.hf(lambda1,lambda0,ARL0,N,h.low,h.up)
H_D1=CUSUM.D.hf(c_theta.D[1],lambda0,ARL0,N,h.low,h.up)
H_D2=CUSUM.D.hf(c_theta.D[2],lambda0,ARL0,N,h.low,h.up)
H_AV1=ACUSUM_VN1.h.f(gamma,lambda0,ARL0,N,h.low,h.up)
H_AV2=ACUSUM_VN2.h.f(gamma,lambda0,ARL0,N,h.low,h.up)
H_AV3=ACUSUM_VN3.h.f(gamma,lambda0,ARL0,N,h.low,h.up)
c(H_s,H_D1,H_D2,H_AV1,H_AV2,H_AV3)
#--------------ARL1.F--------------------------
T_length=length(F_Test.D$Total.Number.of.Tests.Performed)
nt=F_Test.D$Total.Number.of.Tests.Performed[(Phase.t+1):T_length]
Xt=F_Test.D$New.Positives[(Phase.t+1):T_length]

CSUM.1.f=function(h,lambda1,lambda0){
  W=n=x=c()
  plot.stat=0;t=0
  
  while(t<(T_length-Phase.t+1) & plot.stat<h){
    
    t=t+1
    n[t]=nt[t]
    x[t]=Xt[t]
    
    if(t==1){W[t]=max(0,x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))
    } else{W[t]=max(0,W[t-1]+x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))}
    
    plot.stat=W[t]
  }
  return(list("t"=t,"stat"=W))
}
CSD.1.f=function(h,theta.D,lambda0){
  W=n=x=c()
  plot.stat=0;t=0
  
  while(t<(T_length-Phase.t+1) & plot.stat<h){
    
    t=t+1
    n[t]=nt[t]
    x[t]=Xt[t]
    
    lambda_star=lambda0+t*theta.D*sqrt(lambda0)
    
    if(t==1){W[t]=max(0,x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))
    } else{W[t]=max(0,W[t-1]+x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))}
    plot.stat=W[t]
  }
  return(list("t"=t,"stat"=W))
}
ASM.11.f=function(h,gamma,lambda0){
  
  W=n=x=y1=c()
  plot.stat=0;t=0
  
  while(t<(T_length-Phase.t+1) & plot.stat<h){
    
    t=t+1
    n[t]=nt[t]
    x[t]=Xt[t]
    
    if(t==1){
      y1[t]=max(lambda0,(1-gamma)*lambda0+gamma*x[t]/n[t])
      W[t]=max(0,x[t]*log(y1[t]/lambda0)-n[t]*(y1[t]-lambda0))
    } else{
      y1[t]=max(lambda0,(1-gamma)* y1[t-1]+gamma*x[t]/n[t])
      W[t]=max(0,W[t-1]+x[t]*log(y1[t]/lambda0)-n[t]*(y1[t]-lambda0))
    }
    plot.stat=W[t]
  }
  return(list("t"=t,"stat"=W))
}
ASM.21.f=function(theta,h,gamma,lambda0){
  
  W=n=x=y1=y2=c()
  plot.stat=0;t=0
  
  while(t<(T_length-Phase.t+1) & plot.stat<h){
    
    t=t+1
    n[t]=nt[t]
    x[t]=Xt[t]
    
    if(t==1){
      y1[t]=(1-gamma)*lambda0+gamma*x[t]/n[t]
      y2[t]=max(lambda0,y1[t])
      W[t]=max(0,x[t]*log(y2[t]/lambda0)-n[t]*(y2[t]-lambda0))
    } else{
      y1[t]=(1-gamma)*y1[t-1]+gamma*x[t]/n[t]
      y2[t]=max(lambda0,y1[t])
      W[t]=max(0,W[t-1]+x[t]*log(y2[t]/lambda0)-n[t]*(y2[t]-lambda0))
    }
    plot.stat=W[t]
  }
  return(list("t"=t,"stat"=W))
}
ASM.31.f=function(theta,h,gamma,lambda0){
  
  W=n=x=y3=y4=y_xt=y_nt=c()
  plot.stat=0;t=0
  
  while(t<(T_length-Phase.t+1) & plot.stat<h){
    
    t=t+1
    n[t]=nt[t]
    x[t]=Xt[t]
    
    if(t==1){
      y_xt[t]=gamma*x[t]+(1-gamma)*lambda0*n[t]
      y_nt[t]=gamma*n[t]+(1-gamma)*n[t]
      y3[t]=y_xt[t]/y_nt[t]
      y4[t]=max(lambda0,y3[t])
      W[t]=max(0,x[t]*log(y4[t]/lambda0)-n[t]*(y4[t]-lambda0))
    } else{
      y_xt[t]=gamma*x[t]+(1-gamma)*y_xt[t-1]
      y_nt[t]=gamma*n[t]+(1-gamma)*y_nt[t-1]
      y3[t]=y_xt[t]/y_nt[t]
      y4[t]=max(lambda0,y3[t])
      W[t]=max(0,W[t-1]+x[t]*log(y4[t]/lambda0)-n[t]*(y4[t]-lambda0))
    }
    plot.stat=W[t]
  }
  return(list("t"=t,"stat"=W))
}

RL1.S=CSUM.1.f(H_s,lambda1,lambda0)
RL1.D1=CSD.1.f(H_D1,c_theta.D[1],lambda0)
RL1.D2=CSD.1.f(H_D2,c_theta.D[2],lambda0)
RL1.AV1=ASM.11.f(H_AV1,gamma,lambda0)
RL1.AV2=ASM.21.f(theta,H_AV2,gamma,lambda0)
RL1.AV3=ASM.31.f(theta,H_AV3,gamma,lambda0)
c(RL1.S$t,RL1.D1$t,RL1.D2$t,RL1.AV1$t,RL1.AV2$t,RL1.AV3$t)
stopCluster(cl)

#--------------- Control chart plot-------------------
Ws=Wd1=Wd2=Wv1=Wv2=Wv3=n=x=c()
plot.stat=0;t=0
while(t<(T_length-Phase.t+1)){
  
  t=t+1
  n[t]=nt[t]
  x[t]=Xt[t]
  
  if(t==1){Ws[t]=max(0,x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))
  } else{Ws[t]=max(0,Ws[t-1]+x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))}
  plot.stat=Ws[t]
}
n=x==c();plot.stat=0;t=0;theta.D=c_theta.D[1]
while(t<(T_length-Phase.t+1)){
  
  t=t+1
  n[t]=nt[t]
  x[t]=Xt[t]
  
  lambda_star=lambda0+t*theta.D*sqrt(lambda0)
  
  if(t==1){Wd1[t]=max(0,x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))
  } else{Wd1[t]=max(0,Wd1[t-1]+x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))}
  plot.stat=Wd1[t]
}
n=x==c();plot.stat=0;t=0;theta.D=c_theta.D[2]
while(t<(T_length-Phase.t+1)){
  
  t=t+1
  n[t]=nt[t]
  x[t]=Xt[t]
  
  lambda_star=lambda0+t*theta.D*sqrt(lambda0)
  
  if(t==1){Wd2[t]=max(0,x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))
  } else{Wd2[t]=max(0,Wd2[t-1]+x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))}
  plot.stat=Wd2[t]
}
n=x=y1=c()
plot.stat=0;t=0
while(t<(T_length-Phase.t+1)){
  
  t=t+1
  n[t]=nt[t]
  x[t]=Xt[t]
  
  if(t==1){
    y1[t]=max(lambda0,(1-gamma)*lambda0+gamma*x[t]/n[t])
    Wv1[t]=max(0,x[t]*log(y1[t]/lambda0)-n[t]*(y1[t]-lambda0))
  } else{
    y1[t]=max(lambda0,(1-gamma)* y1[t-1]+gamma*x[t]/n[t])
    Wv1[t]=max(0,Wv1[t-1]+x[t]*log(y1[t]/lambda0)-n[t]*(y1[t]-lambda0))
  }
  plot.stat=Wv1[t]
}
n=x=y1=y2=c()
plot.stat=0;t=0
while(t<(T_length-Phase.t+1)){
  
  t=t+1
  n[t]=nt[t]
  x[t]=Xt[t]
  
  if(t==1){
    y1[t]=(1-gamma)*lambda0+gamma*x[t]/n[t]
    y2[t]=max(lambda0,y1[t])
    Wv2[t]=max(0,x[t]*log(y2[t]/lambda0)-n[t]*(y2[t]-lambda0))
  } else{
    y1[t]=(1-gamma)*y1[t-1]+gamma*x[t]/n[t]
    y2[t]=max(lambda0,y1[t])
    Wv2[t]=max(0,Wv2[t-1]+x[t]*log(y2[t]/lambda0)-n[t]*(y2[t]-lambda0))
  }
  plot.stat=Wv2[t]
}
n=x=y3=y4=y_xt=y_nt=c()
plot.stat=0;t=0
while(t<(T_length-Phase.t+1)){
  
  t=t+1
  n[t]=nt[t]
  x[t]=Xt[t]
  
  if(t==1){
    y_xt[t]=gamma*x[t]+(1-gamma)*lambda0*n[t]
    y_nt[t]=gamma*n[t]+(1-gamma)*n[t]
    y3[t]=y_xt[t]/y_nt[t]
    y4[t]=max(lambda0,y3[t])
    Wv3[t]=max(0,x[t]*log(y4[t]/lambda0)-n[t]*(y4[t]-lambda0))
  } else{
    y_xt[t]=gamma*x[t]+(1-gamma)*y_xt[t-1]
    y_nt[t]=gamma*n[t]+(1-gamma)*y_nt[t-1]
    y3[t]=y_xt[t]/y_nt[t]
    y4[t]=max(lambda0,y3[t])
    Wv3[t]=max(0,Wv3[t-1]+x[t]*log(y4[t]/lambda0)-n[t]*(y4[t]-lambda0))
  }
  plot.stat=Wv3[t]
}
Ws
Wd1
Wd2
Wv1
Wv2
Wv3

TL=dim(F_Test.D)[1]
par(mfrow=c(3,2),mar=c(3,3,1,1),mgp=c(2,1,0))
plot(c(1:(length(cT)-Phase.t)),Ws[-length(Ws)],type="b",xlab="Time",
     main="a",ylab="CUSUM statistics",lwd=0.7,cex.main=0.8)
abline(h=H_s,lwd=1, col="red", lty=2)
abline(v=which(Ws>H_s)[1],lwd=1, col="blue", lty=4)

plot(c(1:(length(cT)-Phase.t)),Wd1[-length(Wd1)],type="b",xlab="Time",
     main="b",ylab="CUSUM statistics",lwd=0.7,cex.main=0.8)
abline(h=H_D1,lwd=1, col="red", lty=2)
abline(v=which(Wd1>H_D1)[1],lwd=1, col="blue", lty=4)

plot(c(1:(length(cT)-Phase.t)),Wd2[-length(Wd2)],type="b",xlab="Time",
     main="c",ylab="CUSUM statistics",lwd=0.7,cex.main=0.8)
abline(h=H_D2,lwd=1, col="red", lty=2)
abline(v=which(Wd2>H_D2)[1],lwd=1, col="blue", lty=4)

plot(c(1:(length(cT)-Phase.t)),Wv1[-length(Wv1)],type="b",xlab="Time",
     main="d",ylab="ACUSUM statistics",lwd=0.7,cex.main=0.8)
abline(h=H_AV1,lwd=1, col="red", lty=2)
abline(v=which(Wv1>H_AV1)[1],lwd=1, col="blue", lty=4)

plot(c(1:(length(cT)-Phase.t)),Wv2[-length(Wv2)],type="b",xlab="Time",
     main="e",ylab="ACUSUM statistics",lwd=0.7,cex.main=0.8)
abline(h=H_AV2,lwd=1, col="red", lty=2)
abline(v=which(Wv2>H_AV2)[1],lwd=1, col="blue", lty=4)

plot(c(1:(length(cT)-Phase.t)),Wv3[-length(Wv3)],type="b",xlab="Time",
     main="f",ylab="ACUSUM statistics",lwd=0.7,cex.main=0.8)
abline(h=H_AV3,lwd=1, col="red", lty=2)
abline(v=which(Wv3>H_AV3)[1],lwd=1, col="blue", lty=4)


#gamma=0.002; ARL0=370.3; lambda1=max
#new positive=1~7; Sample=20~150;rate=0.01~0.1;N=20000; 
#H_s=3.41308594;
#H_D1=2.80029297;
#H_D2=3.39843750;
#H_AV1=0.44433594;
#H_AV2=0.19287109;
#H_AV30.05828857
#Result:20 26 21 24 19 11