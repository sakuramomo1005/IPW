#####################################################################################
#                               Make tables
# Read in previous stored RData files and draw summary tables
#####################################################################################


setwd() # set the working path where the RData files are stored

# set empty values

BIAS_ind=c()
BIAS_ex=c()
STD_ind=c()
STD_ex=c()
MCSD_ind=c()
MCSD_ex=c()
Cov_ind=c()
Cov_ex=c()
uncov_ind=c()
uncov_ex=c()

for(k in c(10, 25, 50)){
  
  BIASind=c();BIASex=c();STDind=c();STDex=c();MCind=c();MCex=c();
  unind=c();unex=c();covind=c();covex=c()
  
  for(icc in c(0.05,0.1,0.2)){
    for(iccm in c(0,0.1,0.3,0.5,0.7,0.9)){
      # read in data
      names=paste('wrapk',k,icc,iccm,'.RData',sep='')
      load(names) 
      print(names)
      est_ind=result$est_ind 
      est_ex=result$est_ex
      std_ind=result$std_ind
      std_ex=result$std_ex
      warn_ind=result[[5]]
      warn_ex=result$warn_ex
      misp=result$miss
      maxw=result$maxw
      
      # deal about the non-convergence scenarios
      # non-convergences values are set to NA
      for(i in 1:dim(est_ind)[2]){
        est_ind[,i]=ifelse(warn_ind[,i]==0,est_ind[,i],NA)
        std_ind[,i]=ifelse(warn_ind[,i]==0,std_ind[,i],NA)
        est_ex[,i]=ifelse(warn_ex[,i]==0,est_ex[,i],NA)
        std_ex[,i]=ifelse(warn_ex[,i]==0,std_ex[,i],NA)
      }
      
      # calculate the bias, mcsd, and sd
      estind=apply(est_ind,2,mean,na.rm=TRUE);print(estind)
      estex=apply(est_ex,2,mean,na.rm=TRUE)
      mcsdind=apply(est_ind,2,sd,na.rm=TRUE)[2:dim(est_ind)[2]]
      mcsdex=apply(est_ex,2,sd,na.rm=TRUE)[2:dim(est_ind)[2]]
      stdind=apply(std_ind,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
      stdex=apply(std_ex,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
      mean(misp)/5000
      mean(maxw)
      biasind=abs(estind[2:5]-estind[1])
      biasex=abs(estex[2:5]-estex[1])
      trueind=estind[1]
      trueex=estex[1]
      
      # make a table 
      print(biasind)
      BIASind=rbind(BIASind,c(icc,iccm,biasind))
      BIASex=rbind(BIASex,c(icc,iccm,biasex))
      STDind=rbind(STDind,c(icc,iccm,stdind))
      STDex=rbind(STDex,c(icc,iccm,stdex))
      MCind=rbind(MCind,c(icc,iccm,mcsdind))
      MCex=rbind(MCex,c(icc,iccm,mcsdex))
      unind=rbind(unind,c(icc,iccm,sapply(warn_ind,function(x) sum(x!=0))))
      unex=rbind(unex,c(icc,iccm,sapply(warn_ex,function(x) sum(x!=0))))
      
      # calculate the coverage
      cov_ind=rep(0,dim(est_ind)[2])
      cov_ex=rep(0,dim(est_ind)[2])
      for(i in 1:dim(est_ind)[2]){
        m1=est_ind[,i];s1=std_ind[,i]
        m2=est_ex[,i];s2=std_ex[,i]
        cov_ind[i]=sum(((m1-1.96*s1)<trueind) & ((m1+1.96*s1)>trueind),na.rm=TRUE)
        cov_ex[i]=sum(((m2-1.96*s2)<trueex) & ((m2+1.96*s2)>trueex),na.rm=TRUE)
      }
      
      # make a table
      covind=rbind(covind,c(icc,iccm,cov_ind/1000))
      covex=rbind(covex,c(icc,iccm,cov_ex/1000))
    }
  }
  
  # store the previous resutls
  biasind1=BIASind;biasex1=BIASex;stdind1=STDind;stdex1=STDex;mcind1=MCind;mcex1=MCex;
  Uind1=unind;Uex1=unex;Covind1=covind;Covex1=covex
  
  # get the results of MMI
  BIASind=c();BIASex=c();STDind=c();STDex=c();MCind=c();MCex=c();
  unind=c();unex=c();covind=c();covex=c()
  
  for(icc in c(0.05,0.1,0.2)){
    for(iccm in c(0,0.1,0.3,0.5,0.7,0.9)){
      # read in data
      #load in the true results
      names=paste('wrapk',k,icc,iccm,'.RData',sep='')
      load(names)
      est_ind=result$est_ind
      est_ex=result$est_ex
      estind=apply(est_ind,2,mean,na.rm=TRUE)
      estex=apply(est_ex,2,mean,na.rm=TRUE)
      trueind=estind[1]
      trueex=estex[1]
      
      # read in data, load in the MMI results
      names=paste('wrapk_mmi',k,icc,iccm,'.RData',sep='')
      load(names)
      est_ind=result$est_ind
      est_ex=result$est_ex
      std_ind=result$std_ind
      std_ex=result$std_ex
      warn_ind=result[[5]]
      warn_ex=result$warn_ex
      misp=result$miss
      maxw=result$maxw
      
      # deal with the non-convergence
      est_ind=ifelse(warn_ind==0,est_ind,NA)
      std_ind=ifelse(warn_ind==0,std_ind,NA)
      est_ex=ifelse(warn_ex==0,est_ex,NA)
      std_ex=ifelse(warn_ex==0,std_ex,NA)
      
      # calculate the bias, mcsd, sd
      estind=mean(est_ind,na.rm=TRUE)
      estex=mean(est_ex,na.rm=TRUE)
      mcsdind=sd(est_ind,na.rm=TRUE)
      mcsdex=sd(est_ex,na.rm=TRUE)
      stdind=mean(std_ind,na.rm=TRUE)
      stdex=mean(std_ex,na.rm=TRUE)
      
      biasind=abs(estind-trueind)
      biasex=abs(estex-trueex)
      
      # calculate the coverage
      cov_ind=sum(((est_ind-1.96*std_ind)<trueind) &  ((est_ind+1.96*std_ind)>trueind))/1000
      cov_ex=sum(((est_ex-1.96*std_ex)<trueex) &  ((est_ex+1.96*std_ex)>trueex))/1000
      
      # make a table
      BIASind=rbind(BIASind,c(biasind))
      BIASex=rbind(BIASex,c(biasex))
      STDind=rbind(STDind,c(stdind))
      STDex=rbind(STDex,c(stdex))
      MCind=rbind(MCind,c(mcsdind))
      MCex=rbind(MCex,c(mcsdex))
      unind=rbind(unind,c(sum(warn_ind!=0)))
      unex=rbind(unex,c(sum(warn_ex!=0)))
      covind=rbind(covind,c(cov_ind))
      covex=rbind(covex,c(cov_ex))
    }
  }
  
  # combine the results of IPW and MMI
  BIAS_ind0=round(cbind(biasind1,BIASind),3)
  BIAS_ex0=round(cbind(biasex1,BIASex),3)
  STD_ind0=round(cbind(stdind1,STDind),3)
  STD_ex0=round(cbind(stdex1,STDex),3)
  MCSD_ind0=round(cbind(mcind1,MCind),3)
  MCSD_ex0=round(cbind(mcex1,MCex),3)
  Cov_ind0=round(cbind(Covind1[,c(1,2,4:7)],covind),3)
  Cov_ex0=round(cbind(Covex1[,c(1,2,4:7)],covex),3)
  uncov_ind0=round(cbind(Uind1[,c(1,2,4:7)],unind),3)
  uncov_ex0=round(cbind(Uex1[,c(1,2,4:7)],unex),3)
  
  BIAS_ind=rbind(BIAS_ind,BIAS_ind0)
  BIAS_ex=rbind(BIAS_ex,BIAS_ex0)
  STD_ind=rbind(STD_ind,STD_ind0)
  STD_ex=rbind(STD_ex,STD_ex0)
  MCSD_ind=rbind(MCSD_ind,MCSD_ind0)
  MCSD_ex=rbind(MCSD_ex,MCSD_ex0)
  Cov_ind=rbind(Cov_ind,Cov_ind0)
  Cov_ex=rbind(Cov_ex,Cov_ex0)
  uncov_ind=rbind(uncov_ind,uncov_ind0)
  uncov_ex=rbind(uncov_ex,uncov_ex0)
}

# change the table names
colnames(BIAS_ind)=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
colnames(BIAS_ex)=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
colnames(STD_ind)=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
colnames(STD_ex)=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
colnames(MCSD_ind)=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
colnames(MCSD_ex)=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
colnames(Cov_ind)=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
colnames(Cov_ex)=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
colnames(uncov_ind)=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
colnames(uncov_ex)=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')

k=rep(c(10,25,50),each=18)
BIAS_ind=cbind(k,BIAS_ind)
BIAS_ex=cbind(k,BIAS_ex)
STD_ind=cbind(k,STD_ind)
STD_ex=cbind(k,STD_ex)
MCSD_ind=cbind(k,MCSD_ind)
MCSD_ex=cbind(k,MCSD_ex)
Cov_ind=cbind(k,Cov_ind)
Cov_ex=cbind(k,Cov_ex)
uncov_ind=cbind(k,uncov_ind)
uncov_ex=cbind(k,uncov_ex)

kable(BIAS_ind)
kable(BIAS_ex)
kable(STD_ind)
kable(STD_ex)
kable(MCSD_ind)
kable(MCSD_ex)
kable(Cov_ind)
kable(Cov_ex)
kable(uncov_ind)
kable(uncov_ex)

