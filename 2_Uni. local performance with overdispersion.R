#The Final result
rm(list=ls())

install.packages("predint")
install.packages("foreach")
install.packages("doParallel")

library(predint)
library(foreach)
library(doParallel)

cl=makeCluster(detectCores(logical = F))
registerDoParallel(cl)
###-----------Value definition--------------------
rch=read.csv("T1_uni_h.csv")
ii=1;
N=1000
lambda0=1;ARL0=1/0.0027;h.low=0;h.up=30
c_lambda1=c(1.2,1.3,1.4)
c_d.theta=c(0.001,0.0025)
c_rou=c(0,0.01,0.025,0.05,0.1,0.125,0.15,0.2)
c_gamma=c(0.002,0.005,0.01,0.05)
cACUSUM_VN1.h=cACUSUM_VN2.h=cACUSUM_VN3.h=c()
CUSUM.D.h=CUSUM.S.h=c()
T11=T12=T13=T14=matrix(NA,length(c_rou),3)
TWD=matrix(NA,length(c_rou),length(c_d.theta))
TWS=matrix(NA,length(c_rou),length(c_lambda1))


#---------function definiton--------------------

nt.f=function(ii){
  if(ii==1){runif(1,10,20)}
  else{runif(1,10,15)}
}

#--------------ARL1.find----------------------
ACUSUM_VN1.ARL1.f=function(rou,h,gamma,lambda0,ARL0,ii,N){
  
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  AV1.1=function(rou,gamma,lambda0){
    library(predint)

    W=n=y1=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
      if(rou==0){
        x[t]=rpois(1,lambda0*n[t])
      }else{
        x[t]=rqpois(1,lambda0*n[t],1/(1-rou)^2)$y  #Depending on the R version, determine whether you need to add "$y"
      }
      
      if(t==1){
        y1[t]=max(lambda0,(1-gamma)*lambda0+gamma*x[t]/n[t])
        W[t]=max(0,x[t]*log(y1[t]/lambda0)-n[t]*(y1[t]-lambda0))
      } else{
        y1[t]=max(lambda0,(1-gamma)* y1[t-1]+gamma*x[t]/n[t])
        W[t]=max(0,W[t-1]+x[t]*log(y1[t]/lambda0)-n[t]*(y1[t]-lambda0))
      }
      plot.stat=W[t]
    }
    
    return(t)
  }
  
  RL=foreach(i=1:N,.combine = "c")%dopar% AV1.1(rou,gamma,lambda0)
  ARL1=mean(RL)
  return(ARL1) 
}
ACUSUM_VN2.ARL1.f=function(rou,h,gamma,lambda0,ARL0,ii,N){
  library(predint)
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  AV2.1=function(rou,gamma,lambda0){
    W=n=x=y1=y2=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
      if(rou==0){
        x[t]=rpois(1,lambda0*n[t])
      }else{
        x[t]=rqpois(1,lambda0*n[t],1/(1-rou)^2)$y
      }
      
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
    
    return(t)
  }
  RL=foreach(i=1:N,.combine = "c")%dopar% AV2.1(rou,gamma,lambda0)
  ARL1=mean(RL)
  return(ARL1) 
}
ACUSUM_VN3.ARL1.f=function(rou,h,gamma,lambda0,ARL0,ii,N){
  library(predint)
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  AV3.1=function(rou,gamma,lambda0){
    
    W=n=x=y3=y4=y_xt=y_nt=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
      if(rou==0){
        x[t]=rpois(1,lambda0*n[t])
      }else{
        x[t]=rqpois(1,lambda0*n[t],1/(1-rou)^2)$y
      }
      
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
    return(t)
  }
  RL=foreach(i=1:N,.combine = "c")%dopar% AV3.1(rou,gamma,lambda0)
  ARL1=mean(RL)
  return(ARL1) 
}
CUSUMS.ARL1.f=function(rou,h,lambda1,lambda0,ARL0,ii,N){
  
  library(predint)
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  CS.1=function(rou,lambda1,lambda0){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      t=t+1
      n[t]=nt.f(ii)
      if(rou==0){
        x[t]=rpois(1,lambda0*n[t])
      }else{
        x[t]=rqpois(1,lambda0*n[t],1/(1-rou)^2)$y
      }
      
      if(t==1){W[t]=max(0,x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))
      } else{W[t]=max(0,W[t-1]+x[t]*log(lambda1/lambda0)-n[t]*(lambda1-lambda0))}
      plot.stat=W[t]
    }
    return(t)
  }
  
  RL=foreach(i=1:N,.combine = "c")%dopar% CS.1(rou,lambda1,lambda0)
  ARL1=mean(RL)
  return(ARL1) 

}
CUSUMD.ARL1.f=function(rou,h,theta.D,lambda0,ARL0,ii,N){
  library(predint)
  nt.f=function(ii){
    if(ii==1){runif(1,10,20)}
    else{runif(1,10,15)}
  }
  Ch.1=function(rou,theta.D,lambda0){
    W=n=x=c()
    plot.stat=0;t=0
    
    while(t<1000 & plot.stat<h){
      
      t=t+1
      n[t]=nt.f(ii)
      if(rou==0){
        x[t]=rpois(1,lambda0*n[t])
      }else{
        x[t]=rqpois(1,lambda0*n[t],1/(1-rou)^2)$y
      }
      
      lambda_star=lambda0+t*theta.D*sqrt(lambda0)
      
      if(t==1){W[t]=max(0,x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))
      } else{W[t]=max(0,W[t-1]+x[t]*log(lambda_star/lambda0)-n[t]*(lambda_star-lambda0))}
      plot.stat=W[t]
    }
    return(t)
  }
  
  RL=foreach(i=1:N,.combine = "c")%dopar% Ch.1(rou,theta.D,lambda0)
  ARL1=mean(RL)
  return(ARL1) 
  
}

###-----------get the h with different lambda and different gamma-----------
CUSUM.S.h=rch[1:3,2]
CUSUM.D.h=rch[4:5,2]
vh.matrix=matrix(rch[6:17,2],4,3,byrow = T)
cACUSUM_VN1.h=vh.matrix[,1]
cACUSUM_VN2.h=vh.matrix[,2]
cACUSUM_VN3.h=vh.matrix[,3]

#---T11------
for(i in 1:length(c_rou)){
  gamma=c_gamma[1]
  rou=c_rou[i]
  T11[i,1]=ACUSUM_VN1.ARL1.f(rou,cACUSUM_VN1.h[1],gamma,lambda0,ARL0,ii,N)[1]
  T11[i,2]=ACUSUM_VN2.ARL1.f(rou,cACUSUM_VN2.h[1],gamma,lambda0,ARL0,ii,N)[1]
  T11[i,3]=ACUSUM_VN3.ARL1.f(rou,cACUSUM_VN3.h[1],gamma,lambda0,ARL0,ii,N)[1]
}
#---T12------
for(i in 1:length(c_rou)){
  gamma=c_gamma[2]
  rou=c_rou[i]
  T12[i,1]=ACUSUM_VN1.ARL1.f(rou,cACUSUM_VN1.h[2],gamma,lambda0,ARL0,ii,N)[1]
  T12[i,2]=ACUSUM_VN2.ARL1.f(rou,cACUSUM_VN2.h[2],gamma,lambda0,ARL0,ii,N)[1]
  T12[i,3]=ACUSUM_VN3.ARL1.f(rou,cACUSUM_VN3.h[2],gamma,lambda0,ARL0,ii,N)[1]
}
#---T13------
for(i in 1:length(c_rou)){
  gamma=c_gamma[3]
  rou=c_rou[i]
  T13[i,1]=ACUSUM_VN1.ARL1.f(rou,cACUSUM_VN1.h[3],gamma,lambda0,ARL0,ii,N)[1]
  T13[i,2]=ACUSUM_VN2.ARL1.f(rou,cACUSUM_VN2.h[3],gamma,lambda0,ARL0,ii,N)[1]
  T13[i,3]=ACUSUM_VN3.ARL1.f(rou,cACUSUM_VN3.h[3],gamma,lambda0,ARL0,ii,N)[1]
}
#----T14-----------
for(i in 1:length(c_rou)){
  gamma=c_gamma[4]
  rou=c_rou[i]
  T14[i,1]=ACUSUM_VN1.ARL1.f(rou,cACUSUM_VN1.h[4],gamma,lambda0,ARL0,ii,N)[1]
  T14[i,2]=ACUSUM_VN2.ARL1.f(rou,cACUSUM_VN2.h[4],gamma,lambda0,ARL0,ii,N)[1]
  T14[i,3]=ACUSUM_VN3.ARL1.f(rou,cACUSUM_VN3.h[4],gamma,lambda0,ARL0,ii,N)[1]
}
#-----TWS--------
for(i in 1:length(c_rou)){
  rou=c_rou[i]
  TWS[i,1]=CUSUMS.ARL1.f(rou,CUSUM.S.h[1],c_lambda1[1],lambda0,ARL0,ii,N)
  TWS[i,2]=CUSUMS.ARL1.f(rou,CUSUM.S.h[2],c_lambda1[2],lambda0,ARL0,ii,N)
  TWS[i,3]=CUSUMS.ARL1.f(rou,CUSUM.S.h[3],c_lambda1[3],lambda0,ARL0,ii,N)
}
#-----TWD----------
for(i in 1:length(c_rou)){
  rou=c_rou[i]
  TWD[i,1]=CUSUMD.ARL1.f(rou,CUSUM.D.h[1],c_d.theta[1],lambda0,ARL0,ii,N)
  TWD[i,2]=CUSUMD.ARL1.f(rou,CUSUM.D.h[2],c_d.theta[2],lambda0,ARL0,ii,N)
}
stopCluster(cl)

Table11=cbind(TWS,TWD,T11,T12,T13,T14)
rh=rbind(cACUSUM_VN1.h,cACUSUM_VN2.h,cACUSUM_VN3.h)
crh=c(CUSUM.S.h,CUSUM.D.h,c(rh))
Table1=rbind(crh,Table11)
rownames(Table1)=c("h","0","0.01","0.025","0.05","0.1","0.125","0.15","0.2")
colnames(Table1)=c("CUSUMS","CUSUMS","CUSUMS","CUSUMD","CUSUMD","VN1","VN2","VN3","VN1","VN2","VN3","VN1","VN2","VN3","VN1","VN2","VN3")
Table1
