rm(list=ls())

install.packages("foreach")
install.packages("parallel")
install.packages("doParallel")

library(foreach)
library(parallel)
library(doParallel)

cl=makeCluster(detectCores(logical = F))
registerDoParallel(cl)
###-----------Value definition--------------------
ii=1; 
# if ii=2, it is the code of scenario (3), else if ii=3,it is the code of scenario (4).
N=20000
lambda0=1;ARL0=1/0.0027;h.low=0;h.up=20
c_lambda1=c(1.2,1.3,1.4)
c_d.theta=c(0.001,0.0025)
c_gamma=c(0.002,0.005,0.01,0.05)
ctheta=seq(0,0.5,0.05)
cACUSUM_VN1.h=cACUSUM_VN2.h=cACUSUM_VN3.h=c()
CUSUM.D.h=CUSUM.S.h=c()
TWD=matrix(NA,length(ctheta),length(c_d.theta))
TWS=matrix(NA,length(ctheta),length(c_lambda1))
T11=T12=T13=T14=matrix(NA,length(ctheta),3)

#---------function definiton--------------------
c1=13.8065;c2=11.8532;c3=26.4037
nt.f=function(ii,t){
  if(ii==1){
    nt=c1/(1+exp(-(t-c2)/c3))
  }
  else if(ii==2){
    a=2*c1
    b=1+exp(-(t-(c2+26))/c3)
    nt=a/b
  }else{
    c=c1/2.4
    d=1+exp((t-c2)/c3)
    nt=1+c/d
  }
  return(nt)
}

#---------------h.find--------------------
CUSUM.S.hf=function(lambda1,lambda0,ARL0,ii,N,h.low,h.high){
  c1=13.8065;c2=11.8532;c3=26.4037
  nt.f=function(ii,t){
    if(ii==1){
      nt=c1/(1+exp(-(t-c2)/c3))
    }
    else if(ii==2){
      a=2*c1
      b=1+exp(-(t-(c2+26))/c3)
      nt=a/b
    }else{
      c=c1/2.4
      d=1+exp((t-c2)/c3)
      nt=1+c/d
    }
    return(nt)
  }
  
  CSM.h=function(lambda1,lambda0,ii){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      t=t+1
      n[t]=nt.f(ii,t)
      x[t]=rpois(1,lambda0*n[t])
      if(t==1){W[t]=max(0,x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))
      } else{W[t]=max(0,W[t-1]+x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))}
      plot.stat=W[t]
    }
    t
  }
  
  for(i in 1:100){
    
    h=(h.high+h.low)/2
    RL=foreach(i=1:N,.combine="c")%dopar% CSM.h(lambda1,lambda0,ii)
    
    
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
CUSUM.D.hf=function(theta.D,lambda0,ARL0,ii,N,h.low,h.high){
  c1=13.8065;c2=11.8532;c3=26.4037
  
  nt.f=function(ii,t){
    if(ii==1){
      nt=c1/(1+exp(-(t-c2)/c3))
    }
    else if(ii==2){
      a=2*c1
      b=1+exp(-(t-(c2+26))/c3)
      nt=a/b
    }else{
      c=c1/2.4
      d=1+exp((t-c2)/c3)
      nt=1+c/d
    }
    return(nt)
  }
  
  CSD.h=function(theta.D,lambda0,ii){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      t=t+1
      n[t]=nt.f(ii,t)
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
    RL=foreach(i=1:N,.combine="c")%dopar% CSD.h(theta.D,lambda0,ii)
    
   
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
ACUSUM_VN1.h.f=function(gamma,lambda0,ARL0,ii,N,h.low,h.high){
  c1=13.8065;c2=11.8532;c3=26.4037
  nt.f=function(ii,t){
    if(ii==1){
      nt=c1/(1+exp(-(t-c2)/c3))
    }
    else if(ii==2){
      a=2*c1
      b=1+exp(-(t-(c2+26))/c3)
      nt=a/b
    }else{
      c=c1/2.4
      d=1+exp((t-c2)/c3)
      nt=1+c/d
    }
    return(nt)
  }
  ACSM.1=function(gamma,lambda0,ii){
    
    W=n=x=y1=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii,t)
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
    RL=foreach(i=1:N,.combine="c")%dopar%ACSM.1(gamma,lambda0,ii)
    
   
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
ACUSUM_VN2.h.f=function(gamma,lambda0,ARL0,ii,N,h.low,h.high){
  c1=13.8065;c2=11.8532;c3=26.4037
  nt.f=function(ii,t){
    if(ii==1){
      nt=c1/(1+exp(-(t-c2)/c3))
    }
    else if(ii==2){
      a=2*c1
      b=1+exp(-(t-(c2+26))/c3)
      nt=a/b
    }else{
      c=c1/2.4
      d=1+exp((t-c2)/c3)
      nt=1+c/d
    }
    return(nt)
  }
  ASM.2=function(gamma,lambda0,ii){
    
    W=n=x=y1=y2=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii,t)
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
    RL=foreach(i=1:N,.combine="c")%dopar% ASM.2(gamma,lambda0,ii)
    
    
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
ACUSUM_VN3.h.f=function(gamma,lambda0,ARL0,ii,N,h.low,h.high){
  c1=13.8065;c2=11.8532;c3=26.4037
  nt.f=function(ii,t){
    if(ii==1){
      nt=c1/(1+exp(-(t-c2)/c3))
    }
    else if(ii==2){
      a=2*c1
      b=1+exp(-(t-(c2+26))/c3)
      nt=a/b
    }else{
      c=c1/2.4
      d=1+exp((t-c2)/c3)
      nt=1+c/d
    }
    return(nt)
  }
  ASM.3=function(gamma,lambda0,ii){
    
    W=n=x=y3=y4=y_xt=y_nt=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii,t)
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
    RL=foreach(i=1:N,.combine="c")%dopar%ASM.3(gamma,lambda0,ii) 
    
   
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

#--------------ARL1.find----------------------
CUSUMS.ARL1.f=function(theta,h,lambda1,lambda0,ARL0,ii,N){
  c1=13.8065;c2=11.8532;c3=26.4037
  
  nt.f=function(ii,t){
    if(ii==1){
      nt=c1/(1+exp(-(t-c2)/c3))
    }
    else if(ii==2){
      a=2*c1
      b=1+exp(-(t-(c2+26))/c3)
      nt=a/b
    }else{
      c=c1/2.4
      d=1+exp((t-c2)/c3)
      nt=1+c/d
    }
    return(nt)
  }
  CSUM.1.f=function(theta,h,lambda1,lambda0,ii){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      t=t+1
      n[t]=nt.f(ii,t)
      x[t]=rpois(1,(lambda0+t*theta*sqrt(lambda0))*n[t])
      
      if(t==1){W[t]=max(0,x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))
      } else{W[t]=max(0,W[t-1]+x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))}
      plot.stat=W[t]
    }
    t
  }
  
  RL=foreach(i=1:N,.combine="c")%dopar%CSUM.1.f(theta,h,lambda1,lambda0,ii)
  ARL1=mean(RL)
  return(ARL1)
}
CUSUMD.ARL1.f=function(theta,h,theta.D,lambda0,ARL0,ii,N){
  c1=13.8065;c2=11.8532;c3=26.4037
  nt.f=function(ii,t){
    if(ii==1){
      nt=c1/(1+exp(-(t-c2)/c3))
    }
    else if(ii==2){
      a=2*c1
      b=1+exp(-(t-(c2+26))/c3)
      nt=a/b
    }else{
      c=c1/2.4
      d=1+exp((t-c2)/c3)
      nt=1+c/d
    }
    return(nt)
  }
  CSD.1.f=function(theta,h,theta.D,lambda0,ii){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii,t)
      x[t]=rpois(1,(lambda0+t*theta*sqrt(lambda0))*n[t])
      
      lambda_star=lambda0+t*theta.D*sqrt(lambda0)
      
      if(t==1){W[t]=max(0,x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))
      } else{W[t]=max(0,W[t-1]+x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))}
      plot.stat=W[t]
    }
    t
  }
  
  RL=foreach(i=1:N,.combine="c")%dopar%CSD.1.f(theta,h,theta.D,lambda0,ii)
  ARL1=mean(RL)
  return(ARL1)
}
ACUSUM_VN1.ARL1.f=function(theta,h,gamma,lambda0,ARL0,ii,N){
  c1=13.8065;c2=11.8532;c3=26.4037
  nt.f=function(ii,t){
    if(ii==1){
      nt=c1/(1+exp(-(t-c2)/c3))
    }
    else if(ii==2){
      a=2*c1
      b=1+exp(-(t-(c2+26))/c3)
      nt=a/b
    }else{
      c=c1/2.4
      d=1+exp((t-c2)/c3)
      nt=1+c/d
    }
    return(nt)
  }
  ASM.11.f=function(theta,h,gamma,lambda0,ii){
    
    W=n=x=y1=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii,t)
      x[t]=rpois(1,(lambda0+t*theta*sqrt(lambda0))*n[t])
      
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
  RL=foreach(i=1:N,.combine="c")%dopar%ASM.11.f(theta,h,gamma,lambda0,ii)
  
  
  ARL1=mean(RL)
  return(ARL1) 
}
ACUSUM_VN2.ARL1.f=function(theta,h,gamma,lambda0,ARL0,ii,N){
  c1=13.8065;c2=11.8532;c3=26.4037
   nt.f=function(ii,t){
    if(ii==1){
      nt=c1/(1+exp(-(t-c2)/c3))
    }
    else if(ii==2){
      a=2*c1
      b=1+exp(-(t-(c2+26))/c3)
      nt=a/b
    }else{
      c=c1/2.4
      d=1+exp((t-c2)/c3)
      nt=1+c/d
    }
    return(nt)
  }
   ASM.21.f=function(theta,h,gamma,lambda0,ii){
     
     W=n=x=y1=y2=c()
     plot.stat=0;t=0
     
     while(t<1000 & plot.stat<h){
       
       t=t+1
       n[t]=nt.f(ii,t)
       x[t]=rpois(1,(lambda0+t*theta*sqrt(lambda0))*n[t])
       
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
  
   RL=foreach(i=1:N,.combine="c")%dopar%ASM.21.f(theta,h,gamma,lambda0,ii)
  
  ARL1=mean(RL)
  return(ARL1)
}
ACUSUM_VN3.ARL1.f=function(theta,h,gamma,lambda0,ARL0,ii,N){
  c1=13.8065;c2=11.8532;c3=26.4037
  nt.f=function(ii,t){
    if(ii==1){
      nt=c1/(1+exp(-(t-c2)/c3))
    }
    else if(ii==2){
      a=2*c1
      b=1+exp(-(t-(c2+26))/c3)
      nt=a/b
    }else{
      c=c1/2.4
      d=1+exp((t-c2)/c3)
      nt=1+c/d
    }
    return(nt)
  }
  ASM.31.f=function(theta,h,gamma,lambda0,ii){
    
    W=n=x=y3=y4=y_xt=y_nt=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii,t)
      x[t]=rpois(1,(lambda0+t*theta*sqrt(lambda0))*n[t])
      
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
  RL=foreach(i=1:N,.combine="c")%dopar%ASM.31.f(theta,h,gamma,lambda0,ii)
  ARL1=mean(RL)
  return(ARL1)
}


###-----------get the h by the local performance results-----------
rch=read.csv("T2_Inc_h.csv")
CUSUM.S.h=rch[1:3,2]
CUSUM.D.h=rch[4:5,2]								
cACUSUM_VN1.h=rch[c(6,9,12,15),2]
cACUSUM_VN2.h=rch[c(7,10,13,16),2]
cACUSUM_VN3.h=rch[c(8,11,14,17),2]
##-------------get Table 1 ------------------

#---T11------
for(i in 1:length(ctheta)){
  gamma=c_gamma[1]
  theta=ctheta[i]
  T11[i,1]=ACUSUM_VN1.ARL1.f(theta, cACUSUM_VN1.h[1],gamma,lambda0,ARL0,ii,N)
  T11[i,2]=ACUSUM_VN2.ARL1.f(theta,cACUSUM_VN2.h[1],gamma,lambda0,ARL0,ii,N)
  T11[i,3]=ACUSUM_VN3.ARL1.f(theta,cACUSUM_VN3.h[1],gamma,lambda0,ARL0,ii,N)
}
#---T12------
for(i in 1:length(ctheta)){
  gamma=c_gamma[2]
  theta=ctheta[i]
  T12[i,1]=ACUSUM_VN1.ARL1.f(theta,cACUSUM_VN1.h[2],gamma,lambda0,ARL0,ii,N)
  T12[i,2]=ACUSUM_VN2.ARL1.f(theta,cACUSUM_VN2.h[2],gamma,lambda0,ARL0,ii,N)
  T12[i,3]=ACUSUM_VN3.ARL1.f(theta,cACUSUM_VN3.h[2],gamma,lambda0,ARL0,ii,N)
}
#---T13------
for(i in 1:length(ctheta)){
  gamma=c_gamma[3]
  theta=ctheta[i]
  T13[i,1]=ACUSUM_VN1.ARL1.f(theta,cACUSUM_VN1.h[3],gamma,lambda0,ARL0,ii,N)
  T13[i,2]=ACUSUM_VN2.ARL1.f(theta,cACUSUM_VN2.h[3],gamma,lambda0,ARL0,ii,N)
  T13[i,3]=ACUSUM_VN3.ARL1.f(theta,cACUSUM_VN3.h[3],gamma,lambda0,ARL0,ii,N)
}
#---T14------
for(i in 1:length(ctheta)){
  gamma=c_gamma[4]
  theta=ctheta[i]
  T14[i,1]=ACUSUM_VN1.ARL1.f(theta,cACUSUM_VN1.h[4],gamma,lambda0,ARL0,ii,N)[1]
  T14[i,2]=ACUSUM_VN2.ARL1.f(theta,cACUSUM_VN2.h[4],gamma,lambda0,ARL0,ii,N)[1]
  T14[i,3]=ACUSUM_VN3.ARL1.f(theta,cACUSUM_VN3.h[4],gamma,lambda0,ARL0,ii,N)[1]
}
#----------TWS-------------------
for(i in 1:length(ctheta)){
  theta=ctheta[i]
  TWS[i,1]=CUSUMS.ARL1.f(theta,CUSUM.S.h[1],c_lambda1[1],lambda0,ARL0,ii,N)
  TWS[i,2]=CUSUMS.ARL1.f(theta,CUSUM.S.h[2],c_lambda1[2],lambda0,ARL0,ii,N)
  TWS[i,3]=CUSUMS.ARL1.f(theta,CUSUM.S.h[3],c_lambda1[3],lambda0,ARL0,ii,N)
}
##------TWD----------------
for(i in 1:length(ctheta)){
  theta=ctheta[i]
  TWD[i,1]=CUSUMD.ARL1.f(theta,CUSUM.D.h[1],c_d.theta[1],lambda0,ARL0,ii,N)
  TWD[i,2]=CUSUMD.ARL1.f(theta,CUSUM.D.h[2],c_d.theta[2],lambda0,ARL0,ii,N)
}


Table11=cbind(TWS,TWD,T11,T12,T13,T14)
rh=rbind(cACUSUM_VN1.h,cACUSUM_VN2.h,cACUSUM_VN3.h)
crh=c(CUSUM.S.h,CUSUM.D.h,c(rh))
Table2=rbind(crh,Table11)
colnames(Table2)=c("CUSUMS","CUSUMS","CUSUMS","CUSUMD","CUSUMD","VN1","VN2","VN3","VN1","VN2","VN3","VN1","VN2","VN3","VN1","VN2","VN3")
rownames(Table2)=c("h","0.00","0.05","0.10","0.15","0.20","0.25","0.30","0.35","0.40","0.45","0.50")
Table2
stopCluster(cl)