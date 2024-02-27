#The Final result
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
N=20000
lambda0=1;ARL0=1/0.0027;h.low=0;h.up=30
c_lambda1=c(1.2,1.3,1.4)
c_d.theta=c(0.001,0.0025)
c_gamma=c(0.002,0.005,0.01,0.05)
ctheta=c(0.0000,0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5)
cACUSUM_VN1.h=cACUSUM_VN2.h=cACUSUM_VN3.h=c()
CUSUM.D.h=CUSUM.S.h=c()
TWD=matrix(NA,length(ctheta),length(c_d.theta))
TWS=matrix(NA,length(ctheta),length(c_lambda1))
T11=T12=T13=T14=matrix(NA,length(ctheta),3)

#---------function definiton--------------------

nt.f=function(ii){
  if(ii==1){runif(1,10,20)}
  else{runif(1,10,15)}
}

#---------------h.find--------------------
CUSUM.S.hf=function(lambda1,lambda0,ARL0,ii,N,h.low,h.high){
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  
  CS.h=function(lambda1,lamabda0){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      t=t+1
      n[t]=nt.f(ii)
      x[t]=rpois(1,lambda0*n[t])
      if(t==1){W[t]=max(0,x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))
      } else{W[t]=max(0,W[t-1]+x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))}
      plot.stat=W[t]
    }
    t
  }
  
  for(i in 1:100){
    h=(h.high+h.low)/2
    RL=foreach(i=1:N,.combine="c")%dopar% CS.h(lambda1,lamabda0)
    #ARL=mean(RL)
    
    if(abs(mean(RL)-ARL0)<2 | h.high-h.low<0.00005){
      return(h)
      break
    }else if(mean(RL)-ARL0>2){
      h.high=h
    }else{h.low=h}
    
  }
  return(h)
}
CUSUM.D.hf=function(theta.D,lambda0,ARL0,ii,N,h.low,h.high){
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  
  CD.h=function(theta.D,lambda0){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      t=t+1
      n[t]=nt.f(ii)
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
    RL=foreach(i=1:N,.combine="c")%dopar% CD.h(theta.D,lambda0)
    
    
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
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  ACUM.V1.h=function(gamma,lambda0){
    W=n=x=y1=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
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
    RL=foreach(i=1:N,.combine="c")%dopar%ACUM.V1.h(gamma,lambda0)
    
    
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
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  
  ACUM.V2.h=function(gamma,lambda0){
    
    W=n=x=y1=y2=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
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
    
    RL=foreach(i=1:N,.combine="c")%dopar%ACUM.V2.h(gamma,lambda0)
    
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
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  
  ACUSM.V3.h=function(gamma,lambda0){
    
    W=n=x=y3=y4=y_xt=y_nt=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
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
    RL=foreach(i=1:N,.combine="c")%dopar%ACUSM.V3.h(gamma,lambda0)
    
    
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
ACUSUM_VN1.ARL1.f=function(theta,h,gamma,lambda0,ARL0,ii,N){
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  
  ACU.ARL1.f=function(theta,h,gamma,lambda0,ii){
    W=n=x=y1=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
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
  
  RL=foreach(i=1:N,.combine="c")%dopar%ACU.ARL1.f(theta,h,gamma,lambda0,ii)
  ARL1=mean(RL)
  return(ARL1) 
}
ACUSUM_VN2.ARL1.f=function(theta,h,gamma,lambda0,ARL0,ii,N){
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  ACUM2.ARL1.f=function(theta,h,gamma,lambda0,ii){
    
    W=n=x=y1=y2=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
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
  
  RL=foreach(i=1:N,.combine="c")%dopar%ACUM2.ARL1.f(theta,h,gamma,lambda0,ii)
  
  ARL1=mean(RL)
  return(ARL1) 
}
ACUSUM_VN3.ARL1.f=function(theta,h,gamma,lambda0,ARL0,ii,N){
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  
  ACUM3.ARL1.h=function(theta,h,gamma,lambda0,ii){
    
    W=n=x=y3=y4=y_xt=y_nt=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
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
  RL=foreach(i=1:N,.combine="c")%dopar%ACUM3.ARL1.h(theta,h,gamma,lambda0,ii)
  
  ARL1=mean(RL)
  return(ARL1) 
}
CUSUMS.ARL1.f=function(theta,h,lambda1,lambda0,ARL0,ii,N){
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  
  CS.ARL1.h=function(theta,h,lambda1,lambda0,ii){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      t=t+1
      n[t]=nt.f(ii)
      x[t]=rpois(1,(lambda0+t*theta*sqrt(lambda0))*n[t])
      
      if(t==1){W[t]=max(0,x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))
      } else{W[t]=max(0,W[t-1]+x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))}
      plot.stat=W[t]
    }
    t
  }
  
  RL=foreach(i=1:N,.combine="c")%dopar% CS.ARL1.h(theta,h,lambda1,lambda0,ii)
  ARL1=mean(RL)
  return(ARL1)
}
CUSUMD.ARL1.f=function(theta,h,theta.D,lambda0,ARL0,ii,N){
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  
  CUD.ARL1.f=function(theta,h,theta.D,lambda0,ii){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
      x[t]=rpois(1,(lambda0+t*theta*sqrt(lambda0))*n[t])
      
      lambda_star=lambda0+t*theta.D*sqrt(lambda0)
      
      if(t==1){W[t]=max(0,x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))
      } else{W[t]=max(0,W[t-1]+x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))}
      plot.stat=W[t]
    }
    t
  }
  RL=foreach(i=1:N,.combine="c")%dopar%CUD.ARL1.f(theta,h,theta.D,lambda0,ii)
  
  ARL1=mean(RL)
  return(ARL1)
}


###-----------get the h with different lambda and different gamma-----------
for(i in 1:length(c_gamma)){
  gamma=c_gamma[i]
  cACUSUM_VN1.h[i]=ACUSUM_VN1.h.f(gamma,lambda0,ARL0,ii,N,h.low,h.up)[1]
  cACUSUM_VN2.h[i]=ACUSUM_VN2.h.f(gamma,lambda0,ARL0,ii,N,h.low,h.up)[1]
  cACUSUM_VN3.h[i]=ACUSUM_VN3.h.f(gamma,lambda0,ARL0,ii,N,h.low,h.up)[1]
}
for(i in 1:length(c_lambda1)){
  lambda1=c_lambda1[i]
  CUSUM.S.h[i]=CUSUM.S.hf(lambda1,lambda0,ARL0,ii,N,h.low,h.up)[1]
}
for(i in 1:length(c_d.theta)){
  d.theta=c_d.theta[i]
  CUSUM.D.h[i]=CUSUM.D.hf(d.theta,lambda0,ARL0,ii,N,h.low,h.up)[1]
}

##-------------get Table 1 ------------------

#---T11------
for(i in 1:length(ctheta)){
  gamma=c_gamma[1]
  theta=ctheta[i]
  T11[i,1]=ACUSUM_VN1.ARL1.f(theta, cACUSUM_VN1.h[1],gamma,lambda0,ARL0,ii,N)[1]
  T11[i,2]=ACUSUM_VN2.ARL1.f(theta,cACUSUM_VN2.h[1],gamma,lambda0,ARL0,ii,N)[1]
  T11[i,3]=ACUSUM_VN3.ARL1.f(theta,cACUSUM_VN3.h[1],gamma,lambda0,ARL0,ii,N)[1]
}
#---T12------
for(i in 1:length(ctheta)){
  gamma=c_gamma[2]
  theta=ctheta[i]
  T12[i,1]=ACUSUM_VN1.ARL1.f(theta,cACUSUM_VN1.h[2],gamma,lambda0,ARL0,ii,N)[1]
  T12[i,2]=ACUSUM_VN2.ARL1.f(theta,cACUSUM_VN2.h[2],gamma,lambda0,ARL0,ii,N)[1]
  T12[i,3]=ACUSUM_VN3.ARL1.f(theta,cACUSUM_VN3.h[2],gamma,lambda0,ARL0,ii,N)[1]
}
#---T13------
for(i in 1:length(ctheta)){
  gamma=c_gamma[3]
  theta=ctheta[i]
  T13[i,1]=ACUSUM_VN1.ARL1.f(theta,cACUSUM_VN1.h[3],gamma,lambda0,ARL0,ii,N)[1]
  T13[i,2]=ACUSUM_VN2.ARL1.f(theta,cACUSUM_VN2.h[3],gamma,lambda0,ARL0,ii,N)[1]
  T13[i,3]=ACUSUM_VN3.ARL1.f(theta,cACUSUM_VN3.h[3],gamma,lambda0,ARL0,ii,N)[1]
}
#----T14-----------
for(i in 1:length(ctheta)){
  gamma=c_gamma[4]
  theta=ctheta[i]
  T14[i,1]=ACUSUM_VN1.ARL1.f(theta,cACUSUM_VN1.h[4],gamma,lambda0,ARL0,ii,N)[1]
  T14[i,2]=ACUSUM_VN2.ARL1.f(theta,cACUSUM_VN2.h[4],gamma,lambda0,ARL0,ii,N)[1]
  T14[i,3]=ACUSUM_VN3.ARL1.f(theta,cACUSUM_VN3.h[4],gamma,lambda0,ARL0,ii,N)[1]
}
#-----TWS----------
for(i in 1:length(ctheta)){
  theta=ctheta[i]
  TWS[i,1]=CUSUMS.ARL1.f(theta,CUSUM.S.h[1],c_lambda1[1],lambda0,ARL0,ii,N)
  TWS[i,2]=CUSUMS.ARL1.f(theta,CUSUM.S.h[2],c_lambda1[2],lambda0,ARL0,ii,N)
  TWS[i,3]=CUSUMS.ARL1.f(theta,CUSUM.S.h[3],c_lambda1[3],lambda0,ARL0,ii,N)
}
#------TWD----------
for(i in 1:length(ctheta)){
  theta=ctheta[i]
  TWD[i,1]=CUSUMD.ARL1.f(theta,CUSUM.D.h[1],c_d.theta[1],lambda0,ARL0,ii,N)
  TWD[i,2]=CUSUMD.ARL1.f(theta,CUSUM.D.h[2],c_d.theta[2],lambda0,ARL0,ii,N)
}
stopCluster(cl)

Table11=cbind(TWS,TWD,T11,T12,T13,T14)
rh=rbind(cACUSUM_VN1.h,cACUSUM_VN2.h,cACUSUM_VN3.h)
crh=c(CUSUM.S.h,CUSUM.D.h,c(rh))
Table1=rbind(crh,Table11)

colnames(Table1)=c("CUSUMS","CUSUMS","CUSUMS","CUSUMD","CUSUMD","VN1","VN2","VN3","VN1","VN2","VN3","VN1","VN2","VN3","VN1","VN2","VN3")
rownames(Table1)=c("h","0.0000","0.0001","0.0002","0.0005","0.001","0.002","0.005","0.01","0.02","0.05","0.1","0.2","0.5")
Table1
write.csv(crh,file="T1_uni_h.csv")

