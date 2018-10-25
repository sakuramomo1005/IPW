setwd('/Users/yaolanqiu/Documents/IPW/IPW_final_1011/figures')

### read data
for(read in 1){
  
  ## draw ipw cra tables
  setwd('/Users/yaolanqiu/Documents/IPW/IPW_final_0922/data/geeM')
  
  # set empty values
  BIAS_ind=c();BIAS_ex=c()
  STD_ind=c();STD_ex=c()
  MD_ind=c();MD_ex=c()
  KC_ind=c();KC_ex=c()
  FG_ind=c();FG_ex=c()
  MCSD_ind=c();MCSD_ex=c()
  ROBUST_ind=c();ROBUST_ex=c()
  Cov_ind=c();Cov_ex=c()
  COV_kc_ind=c();COV_kc_ex=c()
  COV_md_ind=c();COV_md_ex=c()
  COV_fg_ind=c();COV_fg_ex=c()
  uncov_ind=c();uncov_ex=c()
  
  for(k in c(10,25,50)){
    
    BIASind=c();BIASex=c();
    STDind=c();STDex=c();Robustind=c();Robustex=c()
    MCind=c();MCex=c();
    MDind=c();MDex=c();KCind=c();KCex=c();FGind=c();FGex=c();
    Cov_md_ind=c();Cov_md_ex=c();Cov_kc_ind=c();Cov_kc_ex=c();
    Cov_fg_ind=c();Cov_fg_ex=c()
    unind=c();unex=c();covind=c();covex=c()
    
    for(icc in c(0.05, 0.1, 0.2)){
      for(iccm in c(0,0.1,0.3,0.5)){
        # read in data
        names=paste('geem',k,icc,iccm,'.RData',sep='')
        load(names) 
        print(names)
        est_ind=result$est_ind 
        est_ex=result$est_ex
        std_ind=result$std_ind
        std_ex=result$std_ex
        
        robust_ind=result$robust_ind
        robust_ex=result$robust_ex
        
        md_ind=result$MD_ind
        md_ex=result$MD_ex
        kc_ind=result$KC_ind
        kc_ex=result$KC_ex
        fg_ind=result$FG_ind
        fg_ex=result$FG_ex
        
        # calculate the bias, mcsd, and sd
        estind=apply(est_ind,2,mean,na.rm=TRUE);print(estind)
        estex=apply(est_ex,2,mean,na.rm=TRUE);print(estex)
        mcsdind=apply(est_ind,2,sd,na.rm=TRUE)[2:dim(est_ind)[2]]
        mcsdex=apply(est_ex,2,sd,na.rm=TRUE)[2:dim(est_ind)[2]]
        stdind=apply(std_ind,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
        stdex=apply(std_ex,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
        robustind=apply(robust_ind,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
        robustex=apply(robust_ex,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
        biasind=100*(estind[2:5]-estind[1])/estind[1]
        biasex=100*(estex[2:5]-estex[1])/estex[1]
        trueind=estind[1]
        trueex=estex[1]
        kcind=apply(kc_ind,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
        kcex=apply(kc_ex,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
        mdind=apply(md_ind,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
        mdex=apply(md_ex,2,mean,na.rm=TRUE) [2:dim(est_ind)[2]]
        fgind=apply(fg_ind,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
        fgex=apply(fg_ex,2,mean,na.rm=TRUE)[2:dim(est_ind)[2]]
        
        # make a table 
        print(biasind)
        BIASind=rbind(BIASind,c(icc,iccm,biasind))
        BIASex=rbind(BIASex,c(icc,iccm,biasex))
        
        STDind=rbind(STDind,c(icc,iccm,stdind))
        Robustind=rbind(Robustind,c(icc,iccm,robustind))
        
        STDex=rbind(STDex,c(icc,iccm,stdex))
        Robustex=rbind(Robustex,c(icc,iccm,robustex))
        
        MCind=rbind(MCind,c(icc,iccm,mcsdind))
        MCex=rbind(MCex,c(icc,iccm,mcsdex))
        
        MDind=rbind(MDind,c(icc,iccm,mdind))
        MDex=rbind(MDex,c(icc,iccm,mdex))
        
        KCind=rbind(KCind,c(icc,iccm,kcind))
        KCex=rbind(KCex,c(icc,iccm,kcex))
        
        FGind=rbind(FGind,c(icc,iccm,fgind))
        FGex=rbind(FGex,c(icc,iccm,fgex))
        
        #corind=rbind(corind,c(icc,iccm,mdind,kcind,fgind))
        #corex=rbind(corex,c(icc,iccm,mdex,kcex,fgex))
        unind=rbind(unind,c(icc,iccm,sapply(est_ind,function(x) sum(is.na(x)))))
        unex=rbind(unex,c(icc,iccm,sapply(est_ex,function(x) sum(is.na(x)))))
        
        # calculate the coverage
        cov_ind=rep(0,dim(est_ind)[2]-1)
        cov_ex=rep(0,dim(est_ind)[2]-1)
        cov_md_ind=rep(0,dim(est_ind)[2]-1)
        cov_md_ex=rep(0,dim(est_ind)[2]-1)
        cov_kc_ind=rep(0,dim(est_ind)[2]-1)
        cov_kc_ex=rep(0,dim(est_ind)[2]-1)
        cov_fg_ind=rep(0,dim(est_ind)[2]-1)
        cov_fg_ex=rep(0,dim(est_ind)[2]-1)
        for(i in 2:dim(est_ind)[2]){
          m1=est_ind[,i];s1=std_ind[,i];s3=md_ind[,i];s5=kc_ind[,i];s7=fg_ind[,i]
          m2=est_ex[,i];s2=std_ex[,i];s4=md_ex[,i];s6=kc_ex[,i];s8=fg_ex[,i]
          cov_ind[i-1]=sum(((m1-1.96*s1)<trueind) & ((m1+1.96*s1)>trueind),na.rm=TRUE)
          cov_ex[i-1]=sum(((m2-1.96*s2)<trueex) & ((m2+1.96*s2)>trueex),na.rm=TRUE)
          
          cov_md_ind[i-1]=sum(((m1-1.96*s3)<trueind) & ((m1+1.96*s3)>trueind),na.rm=TRUE)
          cov_md_ex[i-1]=sum(((m2-1.96*s4)<trueex) & ((m2+1.96*s4)>trueex),na.rm=TRUE)
          
          cov_kc_ind[i-1]=sum(((m1-1.96*s5)<trueind) & ((m1+1.96*s5)>trueind),na.rm=TRUE)
          cov_kc_ex[i-1]=sum(((m2-1.96*s6)<trueex) & ((m2+1.96*s6)>trueex),na.rm=TRUE)
          
          cov_fg_ind[i-1]=sum(((m1-1.96*s7)<trueind) & ((m1+1.96*s7)>trueind),na.rm=TRUE)
          cov_fg_ex[i-1]=sum(((m2-1.96*s8)<trueex) & ((m2+1.96*s8)>trueex),na.rm=TRUE)
        }
        
        # make a table
        covind=rbind(covind,c(icc,iccm,cov_ind/1000))
        covex=rbind(covex,c(icc,iccm,cov_ex/1000))
        
        Cov_md_ind=rbind(Cov_md_ind,c(icc,iccm,cov_md_ind/1000))
        Cov_md_ex=rbind(Cov_md_ex,c(icc,iccm,cov_md_ex/1000))
        
        Cov_kc_ind=rbind(Cov_kc_ind,c(icc,iccm,cov_kc_ind/1000))
        Cov_kc_ex=rbind(Cov_kc_ex,c(icc,iccm,cov_kc_ex/1000))
        
        Cov_fg_ind=rbind(Cov_fg_ind,c(icc,iccm,cov_fg_ind/1000))
        Cov_fg_ex=rbind(Cov_fg_ex,c(icc,iccm,cov_fg_ex/1000))
      }
    }
    
    # change the table names
    colname=c('ICC','ICCM','UCRA','CRA','IPW','IPWC')
    colnames(BIASind)=colnames(BIASex)=
      colnames(STDind)=colnames(STDex)=
      colnames(Robustind)=colnames(Robustex)=
      colnames(KCind)=colnames(KCex)=colnames(MDind)=colnames(MDex)=
      colnames(FGind)=colnames(FGex)=
      colnames(MCind)=colnames(MCex)=
      colnames(Cov_md_ind)=colnames(Cov_md_ex)=
      colnames(Cov_kc_ind)=colnames(Cov_kc_ex)=
      colnames(Cov_fg_ind)=colnames(Cov_fg_ex)=
      colnames(covind)=colnames(covex)=colname
    colnames(unind)=colnames(unex)=c('ICC','ICCM','TRUE','UCRA','CRA','IPW','IPWC')
    
    BIAS_ind=rbind(BIAS_ind,BIASind);BIAS_ex=rbind(BIAS_ex,BIASex);
    MCSD_ind=rbind(MCSD_ind,MCind);MCSD_ex=rbind(MCSD_ex,MCex);
    STD_ind=rbind(STD_ind,STDind);STD_ex=rbind(STD_ex,STDex);
    ROBUST_ind=rbind(ROBUST_ind,Robustind);ROBUST_ex=rbind(ROBUST_ex,Robustex);
    MD_ind=rbind(MD_ind,MDind); MD_ex=rbind(MD_ex,MDex);
    COV_md_ind=rbind(COV_md_ind,Cov_md_ind);COV_md_ex=rbind(COV_md_ex,Cov_md_ex);
    KC_ind=rbind(KC_ind,KCind); KC_ex=rbind(KC_ex,KCex);
    COV_kc_ind=rbind(COV_kc_ind,Cov_kc_ind);COV_kc_ex=rbind(COV_kc_ex,Cov_kc_ex);
    FG_ind=rbind(FG_ind,FGind); FG_ex=rbind(FG_ex,FGex);
    COV_fg_ind=rbind(COV_fg_ind,Cov_fg_ind);COV_fg_ex=rbind(COV_fg_ex,Cov_fg_ex);
    uncov_ind=rbind(uncov_ind,unind);uncov_ex=rbind(uncov_ex,unex);
    Cov_ind=rbind(Cov_ind,covind);Cov_ex=rbind(Cov_ex,covex)
  }
  
  # draw mmi tables
  BIASmmiind=c();BIASmmiex=c()
  STDmmiind=c();STDmmiex=c()
  MDmmiind=c();MDmmiex=c();KCmmiind=c();KCmmiex=c();FGmmiind=c();FGmmiex=c()
  ROBUSTind=c();ROBUSTex=c()
  MCmmiind=c();MCmmiex=c()
  unmmiind=c();unmmiex=c()
  covmmiind=c();covmmiex=c();covMDmmiind=c();covMDmmiex=c()
  covKCmmiind=c();covKCmmiex=c();covFGmmiind=c();covFGmmiex=c()
  
  for(k in c(10,25,50)){
    for(icc in c(0.05,0.1,0.2)){
      for(iccm in c(0,0.1,0.3,0.5)){
        setwd('/Users/yaolanqiu/Documents/IPW/IPW_final_0922/data/geeM')
        
        names=paste('geem',k,icc,iccm,'.RData',sep='')
        load(names) 
        est_ind=result$est_ind 
        est_ex=result$est_ex
        estind=apply(est_ind,2,mean,na.rm=TRUE)
        estex=apply(est_ex,2,mean,na.rm=TRUE)
        trueind=estind[1]
        trueex=estex[1]
        
        names=paste('geem_mmi',k,icc,iccm,'.RData',sep='')
        load(names) 
        print(names)
        est_ind=result$est_ind 
        est_ex=result$est_ex
        std_ind=result$std_ind
        std_ex=result$std_ex
        
        robust_ind=result$ro_ind
        robust_ex=result$ro_ex
        
        md_ind=result$MD_ind
        md_ex=result$MD_ex
        kc_ind=result$KC_ind
        kc_ex=result$KC_ex
        fg_ind=result$FG_ind
        fg_ex=result$FG_ex
        
        df_adj_t = result$df_adj_t
        
        warn_ind=result$warn_ind
        warn_ex=result$warn_ex
        
        # calculate the bias, mcsd, sd
        estind=mean(est_ind,na.rm=TRUE)
        estex=mean(est_ex,na.rm=TRUE)
        mcsdind=sd(est_ind,na.rm=TRUE)
        mcsdex=sd(est_ex,na.rm=TRUE)
        stdind=mean(std_ind,na.rm=TRUE)
        stdex=mean(std_ex,na.rm=TRUE)
        robustind=mean(robust_ind,na.rm=TRUE)
        robustex=mean(robust_ex,na.rm=TRUE)
        
        mdind=mean(md_ind,na.rm=TRUE)
        mdex=mean(md_ex,na.rm=TRUE)
        
        kcind=mean(kc_ind,na.rm=TRUE)
        kcex=mean(kc_ex,na.rm=TRUE)
        
        fgind=mean(fg_ind,na.rm=TRUE)
        fgex=mean(fg_ex,na.rm=TRUE)
        
        biasind=100*(estind-trueind)/trueind
        biasex=100*(estex-trueex)/trueex
        
        # calculate the coverage
        cov_mmi_ind=sum(((est_ind-qt(0.975,df_adj_t)*std_ind)<trueind) &  ((est_ind+qt(0.975,df_adj_t)*std_ind)>trueind))/1000
        cov_mmi_ex=sum(((est_ex-qt(0.975,df_adj_t)*std_ex)<trueex) &  ((est_ex+qt(0.975,df_adj_t)*std_ex)>trueex))/1000
        
        cov_mmi_ro_ind=sum(((est_ind-qt(0.975,df_adj_t)*robust_ind)<trueind) &  ((est_ind+qt(0.975,df_adj_t)*robust_ind)>trueind))/1000
        cov_mmi_ro_ex=sum(((est_ex-qt(0.975,df_adj_t)*robust_ex)<trueex) &  ((est_ex+qt(0.975,df_adj_t)*robust_ex)>trueex))/1000
        
        cov_mmi_MD_ind=sum(((est_ind-qt(0.975,df_adj_t)*md_ind)<trueind) &  ((est_ind+qt(0.975,df_adj_t)*md_ind)>trueind))/1000
        cov_mmi_MD_ex=sum(((est_ex-qt(0.975,df_adj_t)*md_ex)<trueex) &  ((est_ex+qt(0.975,df_adj_t)*md_ex)>trueex))/1000
        
        cov_mmi_KC_ind=sum(((est_ind-qt(0.975,df_adj_t)*kc_ind)<trueind) &  ((est_ind+qt(0.975,df_adj_t)*kc_ind)>trueind))/1000
        cov_mmi_KC_ex=sum(((est_ex-qt(0.975,df_adj_t)*kc_ex)<trueex) &  ((est_ex+qt(0.975,df_adj_t)*kc_ex)>trueex))/1000
        
        cov_mmi_FG_ind=sum(((est_ind-qt(0.975,df_adj_t)*fg_ind)<trueind) &  ((est_ind+qt(0.975,df_adj_t)*fg_ind)>trueind))/1000
        cov_mmi_FG_ex=sum(((est_ex-qt(0.975,df_adj_t)*fg_ex)<trueex) &  ((est_ex+qt(0.975,df_adj_t)*fg_ex)>trueex))/1000
        
        # make a table
        BIASmmiind=rbind(BIASmmiind,c(icc,iccm,biasind))
        BIASmmiex=rbind(BIASmmiex,c(icc,iccm,biasex))
        STDmmiind=rbind(STDmmiind,c(icc,iccm,stdind))
        STDmmiex=rbind(STDmmiex,c(icc,iccm,stdex))
        MDmmiind=rbind(MDmmiind,c(icc,iccm,mdind))
        MDmmiex=rbind(MDmmiex,c(icc,iccm,mdex))
        KCmmiind=rbind(KCmmiind,c(icc,iccm,kcind))
        KCmmiex=rbind(KCmmiex,c(icc,iccm,kcex))
        FGmmiind=rbind(FGmmiind,c(icc,iccm,fgind))
        FGmmiex=rbind(FGmmiex,c(icc,iccm,fgex))
        ROBUSTind=rbind(ROBUSTind,c(icc,iccm,robustind))
        ROBUSTex=rbind(ROBUSTex,c(icc,iccm,robustex))
        MCmmiind=rbind(MCmmiind,c(icc,iccm,mcsdind))
        MCmmiex=rbind(MCmmiex,c(icc,iccm,mcsdex))
        unmmiind=rbind(unmmiind,c(icc,iccm,sum(warn_ind!=0)))
        unmmiex=rbind(unmmiex,c(icc,iccm,sum(warn_ex!=0)))
        covmmiind=rbind(covmmiind,c(icc,iccm,cov_mmi_ind))
        covmmiex=rbind(covmmiex,c(icc,iccm,cov_mmi_ex))
        
        covMDmmiind=rbind(covMDmmiind,c(icc,iccm,cov_mmi_MD_ind))
        covMDmmiex=rbind(covMDmmiex,c(icc,iccm,cov_mmi_MD_ex))
        
        covKCmmiind=rbind(covKCmmiind,c(icc,iccm,cov_mmi_KC_ind))
        covKCmmiex=rbind(covKCmmiex,c(icc,iccm,cov_mmi_KC_ex))
        
        covFGmmiind=rbind(covFGmmiind,c(icc,iccm,cov_mmi_FG_ind))
        covFGmmiex=rbind(covFGmmiex,c(icc,iccm,cov_mmi_FG_ex))
      }
    }
  }
  
  BIASmmiind;BIASmmiex;STDmmiind;STDmmiex;ROBUSTind;MCmmiind;MCmmiex;unmmiind
  unmmiex;covmmiind;covmmiex;covMDmmiind;covMDmmiex;covKCmmiind;covKCmmiex;covFGmmiind;covFGmmiex
  
  
  ## combine ipw, cra and mmi
  BIAS_ind=cbind(BIAS_ind,BIASmmiind[,3])
  BIAS_ex=cbind(BIAS_ex,BIASmmiex[,3])
  STD_ind=cbind(STD_ind,STDmmiind[,3])
  STD_ex=cbind(STD_ex,STDmmiex[,3])
  
  MD_ind=cbind(MD_ind,MDmmiind[,3])
  KC_ind=cbind(KC_ind,KCmmiind[,3])
  FG_ind=cbind(FG_ind,FGmmiind[,3])
  MCSD_ind=cbind(MCSD_ind,MCmmiind[,3])
  ROBUST_ind=cbind(ROBUST_ind,ROBUSTind[,3])
  uncov_ind=cbind(uncov_ind,unmmiind[,3])
  Cov_ind=cbind(Cov_ind,covmmiind[,3])
  COV_md_ind=cbind(COV_md_ind,covMDmmiind[,3])
  COV_kc_ind=cbind(COV_kc_ind,covKCmmiind[,3])
  COV_fg_ind=cbind(COV_fg_ind,covFGmmiind[,3])
  
  MD_ex=cbind(MD_ex,MDmmiex[,3])
  KC_ex=cbind(KC_ex,KCmmiex[,3])
  FG_ex=cbind(FG_ex,FGmmiex[,3])
  MCSD_ex=cbind(MCSD_ex,MCmmiex[,3])
  ROBUST_ex=cbind(ROBUST_ex,ROBUSTex[,3])
  uncov_ex=cbind(uncov_ex,unmmiex[,3])
  Cov_ex=cbind(Cov_ex,covmmiex[,3])
  COV_md_ex=cbind(COV_md_ex,covMDmmiex[,3])
  COV_kc_ex=cbind(COV_kc_ex,covKCmmiex[,3])
  COV_fg_ex=cbind(COV_fg_ex,covFGmmiex[,3])
  
  # change the table names
  colname=c('ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
  colnames(BIAS_ind)=colnames(BIAS_ex)=
    colnames(STD_ind)=colnames(STD_ex)=
    colnames(ROBUST_ind)=colnames(ROBUST_ex)=
    colnames(KC_ind)=colnames(KC_ex)=colnames(MD_ind)=colnames(MD_ex)=
    colnames(FG_ind)=colnames(FG_ex)=
    colnames(MCSD_ind)=colnames(MCSD_ex)=
    colnames(COV_md_ind)=colnames(COV_md_ex)=
    colnames(COV_kc_ind)=colnames(COV_kc_ex)=
    colnames(COV_fg_ind)=colnames(COV_fg_ex)=
    colnames(Cov_ind)=colnames(Cov_ex)=colname
  colnames(uncov_ind)=colnames(uncov_ex)=c('ICC','ICCM','TRUE','UCRA','CRA','IPW','IPWC','MMI')
  
  
  k=c(10,25,50)
  K=rep(k,each=12)
  BIAS_ind=cbind(K,BIAS_ind);BIAS_ex=cbind(K,BIAS_ex)    
  STD_ind=cbind(K,STD_ind);STD_ex=cbind(K,STD_ex) 
  MCSD_ind=cbind(K,MCSD_ind);MCSD_ex=cbind(K,MCSD_ex)
  ROBUST_ind=cbind(K,ROBUST_ind);ROBUST_ex=cbind(K,ROBUST_ex)
  MD_ind=cbind(K,MD_ind);MD_ex=cbind(K,MD_ex)
  KC_ind=cbind(K,KC_ind);KC_ex=cbind(K,KC_ex)
  FG_ind=cbind(K,FG_ind);FG_ex=cbind(K,FG_ex)
  Cov_ind=cbind(K,Cov_ind);Cov_ex=cbind(K,Cov_ex)
  COV_md_ind=cbind(K,COV_md_ind);COV_md_ex=cbind(K,COV_md_ex)
  COV_kc_ind=cbind(K,COV_kc_ind);COV_kc_ex=cbind(K,COV_kc_ex)
  COV_fg_ind=cbind(K,COV_fg_ind);COV_fg_ex=cbind(K,COV_fg_ex)       
  uncov_ind=cbind(K,uncov_ind);uncov_ex=cbind(K,uncov_ex)
  
  
  table_ind=cbind(BIAS_ind,Cov_ind[,4:8],STD_ind[,4:8],MCSD_ind[,4:8],uncov_ind[,5:9])
  table_ex=cbind(BIAS_ex,Cov_ex[,4:8],STD_ex[,4:8],MCSD_ex[,4:8],uncov_ex[,5:9])
}

COV_fg_ind1=COV_fg_ind[K==10,]
COV_md_ind1=COV_md_ind[K==10,]
COV_kc_ind1=COV_kc_ind[K==10,]
Cov_ind1=Cov_ind[K==10,]

COV_fg_ex1=COV_fg_ex[K==10,]
COV_md_ex1=COV_md_ex[K==10,]
COV_kc_ex1=COV_kc_ex[K==10,]
Cov_ex1=Cov_ex[K==10,]

COV=rbind(Cov_ex1,COV_kc_ex1,COV_md_ex1,COV_fg_ex1)
COV=data.frame(COV)
COV$K=rep(c('Cov','KC','MD','FG'),each=12)

Cov_ex[,4:8]=100*Cov_ex[,4:8]
Cov_ind[,4:8]=100*Cov_ind[,4:8]

setwd('/Users/yaolanqiu/Documents/IPW/IPW_final_1011/figures')

# plot figure 2
for(figure in 2){
  # pdf('Figure 2 version 3.pdf')
  jpeg('Figure2.jpeg',width = 7, height = 7, units = 'in', res = 300)
  
  par(oma=c(0,0,2,7))  
  BIAS_ex=data.frame(BIAS_ex)
  colnames(BIAS_ex)=c('K','ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
  ICCM=c(0,0.1,0.3,0.5)
  par(mfrow=c(3,3))
  par(mar=c(4,4,3,1))
  m <- matrix(c(1:9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
  layout(mat=m,heights = c(2,2,2,0.8))
  count=0
  for(icc in c(0.05, 0.1,0.2)){
    for(k in c(10,25,50)){
      count=count+1
      temp=BIAS_ex[BIAS_ex$K==k & BIAS_ex$ICC==icc ,]
      if(count<4){
        kname=paste('K =',k)
      }else{kname=''}
      if(count %%3 ==1){
        iccname='Mean relative bias (%)' #paste('ICC = ',icc)
      }else{iccname=''}
      
      plot(x=ICCM,y=temp$UCRA,type='o',axes = FALSE,
           ylab=iccname,pch=8,cex=1,xlab='ICCM',
           ylim=c(min(BIAS_ex[,4:8]), max(BIAS_ex[,4:8])),main=kname,
           font.main = 2, cex.main = 1.1,lty=1,las=2)
      abline(h=0, col='grey',lwd=1.5,lty = 2)
      lines(ICCM,temp$CRA,col='red',type='o',cex=1,pch=15,lty=5)
      lines(ICCM,temp$IPW,col='blue',type='o',cex=1,pch=16,lty=2)
      lines(ICCM,temp$IPWC,col='darkgreen',type='o',cex=1,pch=17,lty=3)
      lines(ICCM,temp$MMI,col='orange',type='o',cex=1,pch=18,lty=4)
      axis(side=2, at=c(-3,-2,0,2,4,6),las=2)
      axis(side=1, at=c(0,0.1,0.3,0.5))
      if(count %%3==0){
        #  axis(4)
        mtext(paste('ICC = ',icc),font=2,side=4,line = 3,cex=0.75, las=1)
      }
    }
  }
  oldmar<-par(mar=c(1,1,1,1))
  plot.new()
  legend(x = "top",inset = 0,
         legend = c("CRA-GEE", "A-CRA-GEE  ", "W-GEE","CW-GEE", "MMI-GEE"), 
         col=c('black','red','blue','darkgreen','orange'), 
         pch=c(8,15,16,17,18),lty=c(5,1,2,3,4),
         horiz = TRUE)
  par(oldmar)
  dev.off()
}


# plot figure 3
for(figure in 3){
  jpeg('Figure3.jpeg',width = 7, height = 7, units = 'in', res = 300)
  #pdf('Figure 3 version 3.pdf')
  par(oma=c(0,0,2,7))  
  Cov_ex=data.frame(Cov_ex)
  colnames(Cov_ex)=c('K','ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
  ICCM=c(0,0.1,0.3,0.5)
  par(mfrow=c(3,3))
  par(mar=c(4,4,3,1))
  m <- matrix(c(1:9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
  layout(mat=m,heights = c(2,2,2,0.8))
  count=0
  for(icc in c(0.05, 0.1,0.2)){
    for(k in c(10,25,50)){
      count=count+1
      temp=Cov_ex[Cov_ex$K==k & Cov_ex$ICC==icc ,]
      if(count<4){
        kname=paste('K=',k)
      }else{kname=''}
      if(count %%3 ==1){
        iccname='Coverage (%)' #paste('ICC = ',icc)
      }else{iccname=''}
      
      plot(x=ICCM,y=temp$UCRA,type='o',axes = FALSE,
           ylab=iccname,pch=8,cex=1,xlab='ICCM',
           ylim=c(min(Cov_ex[,c(4,6,7,8)]), 100),las=2,
           main=kname,lty=1,font.main=2,cex.main=1.1)
      # lines(ICCM,temp$CRA,col='red',type='o',cex=1,pch=15,lty=5)
      abline(h=95,col='grey',lwd=1.5,lty = 2)
      lines(ICCM,temp$IPW,col='blue',type='o',cex=1,pch=16,lty=2)
      lines(ICCM,temp$IPWC,col='darkgreen',type='o',cex=1,pch=17,lty=3)
      lines(ICCM,temp$MMI,col='orange',type='o',cex=1,pch=18,lty=4)
      axis(side=2, at=c(70,80,85,90,95,100),las=2)
      axis(side=1, at=c(0,0.1,0.3,0.5))
      if(count %%3==0){
        #  axis(4)
        mtext(paste('ICC = ',icc),side=4,line = 3,cex=0.7,font=2,las=1)
      }
    }
  }
  oldmar<-par(mar=c(1,1,1,1))
  plot.new()
  legend(x = "top",inset = 0,
         legend = c("CRA-GEE",  "W-GEE","CW-GEE", "MMI-GEE"), 
         col=c('black','blue','darkgreen','orange'), 
         pch=c(8,16,17,18),lty=c(5,2,3,4),
         horiz = TRUE)
  par(oldmar)
  # title(main="Figure 3: Coverage (%) of four GEE methods to handle missing outcomes in CRTs",outer=T)
  dev.off()
}

# plot figure 4

diff_ex = STD_ex - MCSD_ex
diff_ex[,1:3] = STD_ex[,1:3]

diff_ex = round(diff_ex,3)
for(figure in 4){
  
  jpeg('Figure4.jpeg',width = 7, height = 7, units = 'in', res = 300)
  
  par(oma=c(0,0,2,7))  
  diff_ex=data.frame(diff_ex)
  colnames(diff_ex)=c('K','ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
  ICCM=c(0,0.1,0.3,0.5)
  par(mfrow=c(3,3))
  par(mar=c(4,4,3,1))
  m <- matrix(c(1:9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
  layout(mat=m,heights = c(2,2,2,0.8))
  count=0
  for(icc in c(0.05, 0.1,0.2)){
    for(k in c(10,25,50)){
      count=count+1
      temp=diff_ex[diff_ex$K==k & diff_ex$ICC==icc ,]
      if(count<4){
        kname=paste('K =',k)
      }else{kname=''}
      if(count %%3 ==1){
        iccname='Mean SE - MCSD' #paste('ICC = ',icc)
      }else{iccname=''}
      
      plot(x=ICCM,y=temp$UCRA,type='o',axes = FALSE,
           ylab=iccname,pch=8,cex=1,xlab='ICCM',
           ylim=c(min(diff_ex[,4:8]), max(diff_ex[,4:8])),main=kname,
           font.main = 2, cex.main = 1.1,lty=1,las=2)
      abline(h=0, col='grey',lwd=1.5,lty = 2)
      lines(ICCM,temp$CRA,col='red',type='o',cex=1,pch=15,lty=5)
      lines(ICCM,temp$IPW,col='blue',type='o',cex=1,pch=16,lty=2)
      lines(ICCM,temp$IPWC,col='darkgreen',type='o',cex=1,pch=17,lty=3)
      lines(ICCM,temp$MMI,col='orange',type='o',cex=1,pch=18,lty=4)
      axis(side=2, at=c(-0.3,-0.2,-0.1,0,0.1),las=2)
      axis(side=1, at=c(0,0.1,0.3,0.5))
      if(count %%3==0){
        #  axis(4)
        mtext(paste('ICC = ',icc),font=2,side=4,line = 3,cex=0.75, las=1)
      }
    }
  }
  oldmar<-par(mar=c(1,1,1,1))
  plot.new()
  legend(x = "top",inset = 0,
         legend = c("CRA-GEE", "A-CRA-GEE  ", "W-GEE","CW-GEE", "MMI-GEE"), 
         col=c('black','red','blue','darkgreen','orange'), 
         pch=c(8,15,16,17,18),lty=c(5,1,2,3,4),
         horiz = TRUE)
  par(oldmar)
  dev.off()
}

for(figure in 4){
  
  jpeg('Figure4 without CRA.jpeg',width = 7, height = 7, units = 'in', res = 300)
  
  par(oma=c(0,0,2,7))  
  diff_ex=data.frame(diff_ex)
  colnames(diff_ex)=c('K','ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
  ICCM=c(0,0.1,0.3,0.5)
  par(mfrow=c(3,3))
  par(mar=c(4,4,3,1))
  m <- matrix(c(1:9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
  layout(mat=m,heights = c(2,2,2,0.8))
  count=0
  for(icc in c(0.05, 0.1,0.2)){
    for(k in c(10,25,50)){
      count=count+1
      temp=diff_ex[diff_ex$K==k & diff_ex$ICC==icc ,]
      if(count<4){
        kname=paste('K =',k)
      }else{kname=''}
      if(count %%3 ==1){
        iccname='Mean SE - MCSD' #paste('ICC = ',icc)
      }else{iccname=''}
      
      plot(x=ICCM,y=temp$UCRA,type='o',axes = FALSE,
           ylab=iccname,pch=8,cex=1,xlab='ICCM',
           ylim=c(min(diff_ex[,c(4,6:8)]), max(diff_ex[,c(4,6:8)])),main=kname,
           font.main = 2, cex.main = 1.1,lty=1,las=2)
      abline(h=0, col='grey',lwd=1.5,lty = 2)
     # lines(ICCM,temp$CRA,col='red',type='o',cex=1,pch=15,lty=5)
      lines(ICCM,temp$IPW,col='blue',type='o',cex=1,pch=16,lty=2)
      lines(ICCM,temp$IPWC,col='darkgreen',type='o',cex=1,pch=17,lty=3)
      lines(ICCM,temp$MMI,col='orange',type='o',cex=1,pch=18,lty=4)
      axis(side=2, at=c(-0.2,-0.1,0,0.1),las=2)
      axis(side=1, at=c(0,0.1,0.3,0.5))
      if(count %%3==0){
        #  axis(4)
        mtext(paste('ICC = ',icc),font=2,side=4,line = 3,cex=0.75, las=1)
      }
    }
  }
  oldmar<-par(mar=c(1,1,1,1))
  plot.new()
  legend(x = "top",inset = 0,
         legend = c("CRA-GEE",  "W-GEE","CW-GEE", "MMI-GEE"), 
         col=c('black','blue','darkgreen','orange'), 
         pch=c(8,16,17,18),lty=c(5,2,3,4),
         horiz = TRUE)
  par(oldmar)
  dev.off()
}

# plot figure S1

for(figure in 3){
  jpeg('FigureS1.jpeg',width = 7, height = 7, units = 'in', res = 300)
  #pdf('Figure 3 version 3.pdf')
  par(oma=c(0,0,2,7))  
  Cov_ex=data.frame(Cov_ex)
  colnames(Cov_ex)=c('K','ICC','ICCM','UCRA','CRA','IPW','IPWC','MMI')
  ICCM=c(0,0.1,0.3,0.5)
  par(mfrow=c(3,3))
  par(mar=c(4,4,3,1))
  m <- matrix(c(1:9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
  layout(mat=m,heights = c(2,2,2,0.8))
  count=0
  for(icc in c(0.05, 0.1,0.2)){
    for(k in c(10,25,50)){
      count=count+1
      temp=Cov_ex[Cov_ex$K==k & Cov_ex$ICC==icc ,]
      if(count<4){
        kname=paste('K=',k)
      }else{kname=''}
      if(count %%3 ==1){
        iccname='Coverage (%)' #paste('ICC = ',icc)
      }else{iccname=''}
      
      plot(x=ICCM,y=temp$UCRA,type='o',axes = FALSE,
           ylab=iccname,pch=8,cex=1,xlab='ICCM',
           ylim=c(min(Cov_ex[,c(4,5,6,7,8)]), 100),las=2,
           main=kname,lty=1,font.main=2,cex.main=1.1)
      lines(ICCM,temp$CRA,col='red',type='o',cex=1,pch=15,lty=5)
      abline(h=95,col='grey',lwd=1.5,lty = 2)
      lines(ICCM,temp$IPW,col='blue',type='o',cex=1,pch=16,lty=2)
      lines(ICCM,temp$IPWC,col='darkgreen',type='o',cex=1,pch=17,lty=3)
      lines(ICCM,temp$MMI,col='orange',type='o',cex=1,pch=18,lty=4)
      axis(side=2, at=c(60, 70,80,90,100),las=2)
      axis(side=1, at=c(0,0.1,0.3,0.5))
      if(count %%3==0){
        #  axis(4)
        mtext(paste('ICC = ',icc),side=4,line = 3,cex=0.7,font=2,las=1)
      }
    }
  }
  oldmar<-par(mar=c(1,1,1,1))
  plot.new()
  legend(x = "top",inset = 0,
         legend = c("CRA-GEE", "A-CRA-GEE  ", "W-GEE","CW-GEE", "MMI-GEE"), 
         col=c('black','red','blue','darkgreen','orange'), 
         pch=c(8,15,16,17,18),lty=c(5,1,2,3,4),
         horiz = TRUE)
  par(oldmar)
  # title(main="Figure 3: Coverage (%) of four GEE methods to handle missing outcomes in CRTs",outer=T)
  dev.off()
}
