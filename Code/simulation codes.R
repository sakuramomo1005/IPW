# 2018/08/30

#####################################################################################
#                                    Roadmap
# This is the R code for:
#  1. Generate simulation datasets (datasets with missing outcomes)
#  2. Analyze the datasets with CRA-GEE, A-CRA-GEE,W-GEE,CW-GEE, and MMI-GEE
#. 3. Small sample corrections. 

# Details:
# The simulation steps are:
#  1. Generate a full datasets without missingness
#  2. Generate missing values in outcome based on the function we mention in the paper
#  3. Analyze the datasets; the package for gee is geeM
#  4. Save the results (mean values, sd, mcsd, coverage rate, non-convergence times) in
#     a RData file
#  5. Repeat step 1-4 for 1000 times
#  6. Load in the saved RData files to make a summary table

# MMI
#  1. Impute the dataset for 15 times
#  2. Analyze each imputated dataset with GEE function
#  3. Pool all the results; pool function wroten by myself
#####################################################################################

############### MENU ###############
# 1. Load library
# 2. Functions
  ## 2.1 Function to catch errors and warns
  ## 2.2 Function to calculate the variance based on ICC
  ## 2.3 Function to calculate the missingness percentage in dataset
  ## 2.4 Expit function
  ## 2.5 Data generation function
  ## 2.6 Pool function
  ## 2.7 Small sample corrections (credit to Fan)
# 3. Simulation for MMI-GEE
  ## 3.1. set parameters
  ## 3.2. set empty variables to store the results
  ## 3.3 Tune the intercept
  ## 3.4 Generate a dataset based on the tuned intercept
     ### 3.4.1 Run to generate Nimp full datasets
     ### 3.4.2 Analyze each of the full dataset with GEE
     ### 3.4.3 Pool the results together
# 4. Simulation for CRA-GEE,A-CRA-GEE,W-GEE,and CW-GEE
  ## 4.1. set empty variables to store the results
  ## 4.2 Tune the intercept
  ## 4.3 Generate a dataset based on the tuned intercept
  ## 4.4 Calculate the weights for IPW
  ## 4.5 Calculation
     ### 4.5.1 True effects
     ### 4.5.2 CRA-GEE
     ### 4.5.3 A-CRA-GEE
     ### 4.5.4 W-GEE
     ### 4.5.5 CW-GEE
# 5. Combine results and make tables
####################################

### 1. Load librarys
library(lme4);library(geeM);library(jomo)
library(knitr);library(kableExtra)
setwd("") ## set a dictionary

### 2. Functions
# 2.1 Function to catch errors and warns
myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value, warning=warn, error=err)
}

# 2.2 Function to calculate the variance based on ICC
missing_icc=function(icc){
  pi=3.142
  a=pi^2/(3*(1/icc-1))
  return(a)
}

# 2.3 Function to calculate the missingness percentage in dataset
missing_per=function(data){
  res=sum(data$r)/dim(data)[1]
  return(res)
}

# 2.4 Expit function
expit=function(x){y=exp(x)/(1+exp(x));return(y)}

# 2.5 Data generation function
dategen=function(k,M,mux=0,varx,icc,mud=0,iccm,intercept){
  ## intervention group
  K=2*k  # total cluster number
  m=rpois(K,M) # cluster sizes
  N=sum(m) # total individual number
  i=rep(rep(c(0,1),each=k),times=m) # 
  cluster=rep(1:K,times=m)
  vard=missing_icc(icc)
  delta=rep(rnorm(K,mud,sqrt(vard)),times=m)
  x=rnorm(N,mux,sqrt(varx))
  p=expit(1+1.36*i+x+delta)
  y=rbinom(N,1,p)
  alpha=rep(rnorm(K,0,sqrt(missing_icc(iccm))),times=m)
  mis=expit(intercept+i+x+alpha)
  r=rbinom(N,1,mis)
  res=data.frame(y=y,arm=i,x=x,cluster=cluster,delta=delta,mis=mis,r=r)
  return(res)
}

# 2.6 Pool function
mypool=function(mean0,sd0,num=5,print='no'){
  m=mean(mean0,na.rm=TRUE)
  v=mean(sd0,na.rm=TRUE)
  B=sd(mean0,na.rm=TRUE)
  v_hat=v+(1+1/num)*B
  l=m-1.96*v_hat
  u=m+1.96*v_hat
  if(print=='no'){
    return(list(mean=m,std=v_hat))
  }
  if(print=='yes'){
    print('mean (95% CI)')
    print(paste(round(m,2)," (",round(l,2),',',round(u,2),')',sep=''))
    return(list(mean=m,std=v_hat))
  }
}

# 2.7 Small sample corrections (Fan)
BCVexch=function(y,X,beta,alpha,phi,id,w){
  
  require(MASS)
  
  # Creates two vectors that have the start and end points for each cluster
  BEGINEND=function(n){
    last=cumsum(n)
    first=last-n+1
    return(cbind(first,last))
  }
  
  # Score function
  SCORE=function(beta,alpha,phi,y,X,n,p){
    U=rep(0,p)
    UUtran=Ustar=matrix(0,p,p)
    locx=BEGINEND(n)
    
    for(i in 1:length(n)){
      X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
      y_c=y[locx[i,1]:locx[i,2]]
      w_c=w[locx[i,1]:locx[i,2]]
      
      U_c=rep(0,p)
      Ustar_c=matrix(0,p,p)
      mu_c=1/(1+exp(c(-X_c%*%beta)))
      
      C=X_c*(mu_c*(1-mu_c))
      A=y_c-mu_c
      D=diag(w_c,nrow=length(w_c))
      INVR=diag(1/(1-alpha),n[i])-matrix(alpha/((1-alpha)*(1-alpha+n[i]*alpha)),n[i],n[i])
      INVB=diag(1/sqrt(mu_c*(1-mu_c)),n[i]) %*% INVR %*% diag(1/sqrt(mu_c*(1-mu_c)),n[i]) %*% D/phi
      
      U_c=t(C)%*%INVB%*%A
      UUtran_c=tcrossprod(U_c)
      Ustar_c=t(C)%*%INVB%*%C
      U=U+U_c
      UUtran=UUtran+UUtran_c
      Ustar=Ustar+Ustar_c
    }
    return(list(U=U,UUtran=UUtran,Ustar=Ustar))
  }
  
  # Creates bias-corrected covariance matrix of beta
  p=ncol(X)
  n=as.numeric(table(id))
  SCORE_RES=SCORE(beta,alpha,phi,y,X,n,p)
  U=SCORE_RES$U
  UUtran=SCORE_RES$UUtran
  Ustar=SCORE_RES$Ustar
  
  # Naive or Model-based estimator
  naive=ginv(Ustar)
  
  # BC0 or usual Sandwich estimator     
  robust=naive%*%UUtran%*%t(naive)
  
  # Bias-corrected variance
  Ustar_c_array=UUtran_c_array=array(0,c(p,p,length(n)))
  UUtran=UUbc=UUbc2=UUbc3=Ustar=matrix(0,p,p)
  
  locx=BEGINEND(n)
  
  for(i in 1:length(n)){
    X_c=X[locx[i,1]:locx[i,2],,drop=FALSE]
    y_c=y[locx[i,1]:locx[i,2]]
    w_c=w[locx[i,1]:locx[i,2]]
    mu_c=1/(1+exp(c(-X_c%*%beta)))
    
    U_i=U_c=rep(0,p)
    Ustar_c=matrix(0,p,p)
    
    # commands for beta
    C=X_c*(mu_c*(1-mu_c))
    A=y_c-mu_c
    D=diag(w_c,nrow=length(w_c))
    INVR=diag(1/(1-alpha),n[i])-matrix(alpha/((1-alpha)*(1-alpha+n[i]*alpha)),n[i],n[i])
    INVB=diag(1/sqrt(mu_c*(1-mu_c)),n[i]) %*% INVR %*% diag(1/sqrt(mu_c*(1-mu_c)),n[i]) %*% D/phi
    U_i=t(C)%*%INVB%*%A
    
    L_i=ginv(diag(1,nrow=n[i]) - C%*%naive%*%t(C)%*%INVB) %*% A
    U_c=t(C)%*%INVB%*%L_i
    
    Ustar_c=t(C)%*%INVB%*%D%*%C
    Ustar=Ustar+Ustar_c
    UUtran_c=tcrossprod(U_i)
    UUtran=UUtran+UUtran_c
    UUbc_c=tcrossprod(U_c)
    UUbc=UUbc+UUbc_c
    UUbc_ic=tcrossprod(U_c,U_i)
    UUbc2=UUbc2+UUbc_ic
    
    Ustar_c_array[,,i]=Ustar_c
    UUtran_c_array[,,i]=UUtran_c
  }
  
  # calculating adjustment factor for BC3
  for(i in 1:length(n)){      
    Hi=diag(1/sqrt(1-pmin(0.75,c(diag(Ustar_c_array[,,i]%*%naive)))))
    UUbc3=UUbc3+Hi%*%UUtran_c_array[,,i]%*%Hi
  }
  
  # BC1 or Variance estimator due to Kauermann and Carroll (2001);
  varKC=naive%*%(UUbc2+t(UUbc2))%*%t(naive)/2
  
  # BC2 or Variance estimator due to Mancl and DeRouen (2001);
  varMD=naive%*%UUbc%*%t(naive)
  
  # BC3 or Variance estimator due to Fay and Graubard (2001);
  varFG=naive%*%UUbc3%*%t(naive)
  
  ########################################
  # Output
  # naive: naive or model-based var
  # robust: robust sandwich var
  # varMD: bias-corrected sandwich var due to Mancl and DeRouen (2001)
  # varKC: bias-corrected sandwich var due to Kauermann and Carroll (2001)
  # varFG: bias-corrected sandwich var due to Fay and Graubard (2001)
  ########################################
  return(list(naive=naive,robust=robust,varMD=varMD,varKC=varKC,varFG=varFG))
}

### 3. Simulation for MMI-GEE

# 3.1. set parameters;
varx=0.2; M=50; Nimp=15;
s=Sys.time() # save time 
T=1000 # times of simulation

for(k in c(10, 25, 50)){ # number of clusters 
  for(icc in c(0.05, 0.1, 0.2)){ # icc for data generation
    for(iccm in c(0, 0.1, 0.3, 0.5)){ # icc for missingness
      
      # 3.2. set empty variables to store the results
      est_ind=c();est_ex=c() # Estimates for independent and exchangeable
      std_ind=c();std_ex=c() # SEs for independent and exchangeable
      ro_ind=c();ro_ex=c() # the robust SEs calculated by fan's function, which is same with SEs calculated by geeM
      warn_ind=c();warn_ex=c() # whether there is warning
      MD_ind=c();MD_ex=c();KC_ind=c();KC_ex=c();FG_ind=c();FG_ex=c() # the SEs from each small sample correction methods 
      
      for(times in 1:T){
        print(paste('k',k,'icc',icc,'iccm',iccm,'times',times))
        set.seed(times)
        
        # 3.3 Tune the intercept
        # since the parameters change may bring different missing percentages
        # we tune the intercept to keep the missing percentage around 30%
        
        mis=c()
        for(intercept in seq(-5,5,0.1)){
          temp=dategen(k,M,varx=varx,icc=icc,iccm=iccm,intercept=intercept) 
          mis=c(mis,missing_per(temp))
        }
        intercept=seq(-5,5,0.1)[which.min(abs(mis-0.3))]
        
        # 3.4 generate a dataset based on the tuned intercept
        # d1: the full dataset
        # d3: the dataset with misisng outcomes
        
        d1=dategen(k,M,varx=varx,icc=icc,iccm=iccm,intercept=intercept)
        d3=d1
        d3$y=ifelse(d3$r==1,NA,d3$y)
        d3$missing=d3$r
        
        # 3.4 MMI; just copy Hossain's code
        data.miss=d3
        y.cat= data.frame(outcome=data.miss$y)  # data frame for response variables with missing values
        y.numcat=c(2)                                 # number of levels in outcome variable
        clus=data.frame(clus=data.miss$cluster)          # data frame for clusters
        
        nobs=dim(data.miss)[1]
        x= data.frame(intercept=rep(1,nobs),covariate=data.miss$x,group=data.miss$arm)

        # 3.4.1 run to generate Nimp full datasets
        # "imp" represents full datasets 
        imp = jomo1rancat(Y.cat=y.cat, Y.numcat=y.numcat, X=x,
                          clus=clus,nburn=100, nbetween=25, nimp=Nimp,output=0)
        
        # set empty results and save them later
        mmi_est_ind=c();mmi_std_ind=c();mmi_warn_ind=c()
        mmi_est_ex=c();mmi_std_ex=c();mmi_warn_ex=c()
        mmi_MD_ex=c();mmi_KC_ex=c();mmi_FG_ex=c()
        mmi_FG_ind=c();mmi_KC_ind=c(); mmi_MD_ind=c()
        mmi_ro_ex=c();  mmi_ro_ind=c()
        
        # 3.4.2 Analyze each of the full dataset with GEE
        for(i in 1:Nimp){
          temp=imp[imp$Imputation==i,]
          rownames(temp)=NULL
          temp$outcome=as.numeric(temp$outcome)-1
          
          mmi_ind=myTryCatch(geem(formula=outcome~group,
                                  id=clus , data = temp,
                                  family =  binomial("logit"),
                                  corstr = "independence"))
          
          mmi_ex=myTryCatch(geem(formula=outcome~group,
                                 id=clus , data = temp,
                                 family =  binomial("logit"),
                                 corstr = "exchangeable"))
         
          ## Save the results
          if(is.null(mmi_ind$value)==0){
            phi1=mmi_ind$value$phi
            alpha1=mmi_ind$value$alpha
            beta1=mmi_ind$value$beta
            t1=beta1[2]
            y1=temp$outcome
            X1=cbind(rep(1,length(temp$outcome)),temp$group)
            w1=rep(1,length(temp$outcome))
            id1=temp$clus
            correction1=myTryCatch(BCVexch(y1,X1,beta1,alpha1,phi1,id1,w1))
            std1=summary(mmi_ind$value)$se.robust[2]
            if(is.null(correction1$error)==1){
              correction1=correction1$value
              mmi_est_ind=c(mmi_est_ind,t1)
              mmi_std_ind=c(mmi_std_ind,std1)
              mmi_ro_ind=c(mmi_ro_ind,sqrt(diag(correction1$robust))[2])
              mmi_MD_ind=c(mmi_MD_ind,sqrt(diag(correction1$varMD))[2])
              mmi_KC_ind=c(mmi_KC_ind,sqrt(diag(correction1$varKC))[2])
              mmi_FG_ind=c(mmi_FG_ind,sqrt(diag(correction1$varFG))[2])
              if(is.null(correction1$error)==0){mmi_warn_ind=c(mmi_warn_ind,times)}
              if(is.null(correction1$error)==1){mmi_warn_ind=c(mmi_warn_ind,0)}
            }else{
              mmi_est_ind=c(mmi_est_ind,NA)
              mmi_std_ind=c(mmi_std_ind,NA)
              mmi_ro_ind=c(mmi_ro_ind,NA)
              mmi_MD_ind=c(mmi_MD_ind,NA)
              mmi_KC_ind=c(mmi_KC_ind,NA)
              mmi_FG_ind=c(mmi_FG_ind,NA)
              mmi_warn_ind=c(mmi_warn_ind,times)
            }
          }
          
          if(is.null(mmi_ind$value)==1){
            mmi_est_ind=c(mmi_est_ind,NA)
            mmi_std_ind=c(mmi_std_ind,NA)
            mmi_ro_ind=c(mmi_ro_ind,NA)
            mmi_MD_ind=c(mmi_MD_ind,NA)
            mmi_KC_ind=c(mmi_KC_ind,NA)
            mmi_FG_ind=c(mmi_FG_ind,NA)
            mmi_warn_ind=c(mmi_warn_ind,times)
          }
          
          if( is.null(mmi_ex$value)==0){
            phi2=mmi_ex$value$phi
            alpha2=mmi_ex$value$alpha
            beta2=mmi_ex$value$beta
            t2=beta2[2]
            y2=temp$outcome
            X2=cbind(rep(1,length(temp$outcome)),temp$group)
            w2=rep(1,length(temp$outcome))
            id2=temp$clus
            correction2=myTryCatch(BCVexch(y2,X2,beta2,alpha2,phi2,id2,w2))
            std2=summary(mmi_ex$value)$se.robust[2]
            if(is.null(correction2$error)==1){
              correction2=correction2$value
              mmi_est_ex=c(mmi_est_ind,t2)
              mmi_std_ex=c(mmi_std_ind,std2)
              mmi_ro_ex=c(mmi_ro_ex,sqrt(diag(correction2$robust))[2])
              mmi_MD_ex=c(mmi_MD_ex,sqrt(diag(correction2$varMD))[2])
              mmi_KC_ex=c(mmi_KC_ex,sqrt(diag(correction2$varKC))[2])
              mmi_FG_ex=c(mmi_FG_ex,sqrt(diag(correction2$varFG))[2])
              
              if(is.null(mmi_ex$value)==1){mmi_warn_ex=c(mmi_warn_ex,times)}
              if(is.null(mmi_ex$value)==0){mmi_warn_ex=c(mmi_warn_ex,0)}
            }else{
              mmi_est_ex=c(mmi_est_ex,NA)
              mmi_std_ex=c(mmi_std_ex,NA)
              mmi_ro_ex=c(mmi_ro_ex,NA)
              mmi_MD_ex=c(mmi_MD_ex,NA)
              mmi_KC_ex=c(mmi_KC_ex,NA)
              mmi_FG_ex=c(mmi_FG_ex,NA)
              mmi_warn_ex=c(mmi_warn_ex,times)
            }
          }
          if(is.null(mmi_ex$value)==1){
            mmi_est_ex=c(mmi_est_ex,NA)
            mmi_ro_ex=c(mmi_ro_ex,NA)
            mmi_std_ex=c(mmi_std_ex,NA)
            mmi_MD_ex=c(mmi_MD_ex,NA)
            mmi_KC_ex=c(mmi_KC_ex,NA)
            mmi_FG_ex=c(mmi_FG_ex,NA)
            mmi_warn_ex=c(mmi_warn_ex,times)
          }
        }
        
        # 3.4.3 Pool the results together
        # pool with orginal SE
        temp1=mypool(mmi_est_ind,mmi_std_ind,num=Nimp)
        temp2=mypool(mmi_est_ex,mmi_std_ex,num=Nimp)
        
        temp11=mypool(mmi_est_ind,mmi_ro_ind,num=Nimp)
        temp22=mypool(mmi_est_ex,mmi_ro_ex,num=Nimp)
        
        # pool with corrected SE
        temp3=mypool(mmi_est_ind,mmi_MD_ind,num=Nimp)
        temp4=mypool(mmi_est_ex,mmi_MD_ex,num=Nimp)
        
        temp5=mypool(mmi_est_ind,mmi_KC_ind,num=Nimp)
        temp6=mypool(mmi_est_ex,mmi_KC_ex,num=Nimp)
        
        temp7=mypool(mmi_est_ind,mmi_FG_ind,num=Nimp)
        temp8=mypool(mmi_est_ex,mmi_FG_ex,num=Nimp)
        
        # save results
        est_ind=c(est_ind,temp1$mean)
        std_ind=c(std_ind,temp1$std)
        ro_ind=c(ro_ind,temp11$std)
        ro_ex=c(ro_ex,temp22$std)
        est_ex=c(est_ex,temp2$mean)
        std_ex=c(std_ex,temp2$std)
        
        MD_ind=c(MD_ind,temp3$std)
        MD_ex=c(MD_ex,temp4$std)
        KC_ind=c(KC_ind,temp5$std)
        KC_ex=c(KC_ex,temp6$std)
        FG_ind=c(FG_ind,temp7$std)
        FG_ex=c(FG_ex,temp8$std)
        
        warn_ind=c(warn_ind,sum(mmi_warn_ind))
        warn_ex=c(warn_ex,sum(mmi_warn_ex))
      }
      
      result=list(est_ind=est_ind,est_ex=est_ex,std_ind=std_ind,
                  std_ex=std_ex,ro_ind=ro_ind,
                  ro_ex=ro_ex,warn_ind=warn_ind,warn_ex=warn_ex,
                  MD_ind=MD_ind,MD_ex=MD_ex,KC_ind=KC_ind,KC_ex=KC_ex,
                  FG_ind=FG_ind,FG_ex=FG_ex)
      print(result)
      names=paste('geem_mmi',k,icc,iccm,'.RData',sep='')
      save(result,file=names)
    }
  }
}
e=Sys.time()

### 4. Simulation for CRA, A-CRA,W-GEE,CW-GEE

for(k in c(10, 25, 50)){
  for(icc in c(0.05, 0.1, 0.2)){
    for(iccm in c(0, 0.1, 0.3, 0.5)){
      
      # 4.1 set empty results and save for later
      # the results for true values.
      robust_true_ind=c();MD_true_ind=c();KC_true_ind=c();FG_true_ind=c()
      est_true_ind=c();std_true_ind=c()
      robust_true_ex=c();MD_true_ex=c();KC_true_ex=c();FG_true_ex=c()
      est_true_ex=c();std_true_ex=c()
      
      # the results for CRA-GEE
      robust_ucra_ind=c();MD_ucra_ind=c();KC_ucra_ind=c();FG_ucra_ind=c()
      est_ucra_ind=c();std_ucra_ind=c()
      robust_ucra_ex=c();MD_ucra_ex=c();KC_ucra_ex=c();FG_ucra_ex=c()
      est_ucra_ex=c();std_ucra_ex=c()
      
      # the results for A-CRA-GEE
      robust_cra_ind=c();MD_cra_ind=c();KC_cra_ind=c();FG_cra_ind=c()
      est_cra_ind=c();std_cra_ind=c()
      robust_cra_ex=c();MD_cra_ex=c();KC_cra_ex=c();FG_cra_ex=c()
      est_cra_ex=c();std_cra_ex=c()
      
      # the results for W-GEE
      robust_ipw_ex=c();MD_ipw_ex=c();KC_ipw_ex=c();FG_ipw_ex=c()
      est_ipw_ex=c();std_ipw_ex=c()
      robust_ipw_clu_ind=c();MD_ipw_clu_ind=c();KC_ipw_clu_ind=c();FG_ipw_clu_ind=c()
      est_ipw_clu_ind=c();std_ipw_clu_ind=c()
      
      # the results for CW-GEE
      robust_ipw_clu_ex=c();MD_ipw_clu_ex=c();KC_ipw_clu_ex=c();FG_ipw_clu_ex=c()
      est_ipw_clu_ex=c();std_ipw_clu_ex=c()
      robust_ipw_ind=c();MD_ipw_ind=c();KC_ipw_ind=c();FG_ipw_ind=c()
      est_ipw_ind=c();std_ipw_ind=c()
      
      for(times in 1:1000){
        print(paste('k',k,'icc',icc,'iccm',iccm,'times',times))
        set.seed(times)
        
        ## 4.2 Tune the intercept
        # since the parameters change may bring different missing percentages,
        # we tune the intercept to keep the missing percentage around 30%
        mis=c()
        for(intercept in seq(-5,5,0.1)){
          temp=dategen(k,M,varx=varx,icc=icc,iccm=iccm,intercept=intercept)  
          mis=c(mis,missing_per(temp))
        }
        intercept=seq(-5,5,0.1)[which.min(abs(mis-0.3))]
        
        # 4.3 Generate a dataset based on the tuned intercept
        # d1: the full dataset
        # d3: the dataset with misisng outcomes
        # d2: the d3 dataset that omits the missing values 
        
        d1=dategen(k,M,varx=varx,icc=icc,iccm=iccm,intercept=intercept)
        d3=d1
        d3$y=ifelse(d3$r==1,NA,d3$y)
        d3$missing=d3$r
        
        # 4.4 calculate the weights for IPW
        # w1: the weights without consideration of cluster effects
        # w2: the weights with consideration of cluster effects
        
        w1=glm(missing ~x + arm, data = d3,
               family = binomial(link='logit'))
        w2=glmer(missing ~x + arm+(1|cluster) , data = d3,
                 family = binomial(link='logit'))

        w1=predict(w1,type="response")  # get the weights value from the glm
        w2=predict(w2,type="response")  # get the weights value form the glmer
        
        d3$weight=1/w1
        d3$weight2=1/w2
        
        d2=na.omit(d3)
        
        # 4.5 Calculation
        ### 4.5.1 True effect
        trues_ind=myTryCatch(geem(formula=y~arm,id=cluster, data = d1,
                                  family =  binomial("logit"),
                                  corstr = "independence"))
        trues_ex=myTryCatch(geem(formula=y~arm,id=cluster, data = d1,
                                 family =  binomial("logit"),
                                 corstr = "exchangeable"))
        
        ### 4.5.2 Unadjusted CRA
        ucra_ind=myTryCatch(geem(formula=y~arm,id=cluster, data = d2,
                                 family =  binomial("logit"),
                                 corstr = "independence"))
        
        ucra_ex=myTryCatch(geem(formula=y~arm,id=cluster, data = d2,
                                family =  binomial("logit"),
                                corstr = "exchangeable"))
        
        ### 4.5.3 Adjusted CRA
        cra_ind=myTryCatch(geem(formula=y~x+arm,id=cluster, data = d2,
                                family =  binomial("logit"),
                                corstr = "independence"))
        cra_ex=myTryCatch(geem(formula=y~x+arm,id=cluster, data = d2,
                               family =  binomial("logit"),
                               corstr = "exchangeable"))
        
        ### 4.5.4 IPW without cluster effects 
        ipw_ind=myTryCatch(geem(formula=y~arm,id=cluster, data = d3,
                                family =  binomial("logit"),
                                weights = d3$weight,
                                corstr = "independence"))
        ipw_ex=myTryCatch(geem(formula=y~arm,id=cluster, data = d3,
                               family =  binomial("logit"),
                               weights = d3$weight,
                               corstr = "exchangeable"))
        
        ### 4.5.5 IPW with cluster effects
        ipw_clu_ind=myTryCatch(geem(formula=y~arm,id=cluster, data = d4,
                                    family =  binomial("logit"),
                                    weights = d4$weight2,
                                    corstr = "independence"))
        ipw_clu_ex=myTryCatch(geem(formula=y~arm,id=cluster, data = d3,
                                   family =  binomial("logit"),
                                   weights = d3$weight2,
                                   corstr = "exchangeable"))
        
        # 4.6 Save the results
        ### true independent
        if(is.null(trues_ind$value)==0){
          if(trues_ind$value$converged==1){
            phi1=trues_ind$value$phi
            alpha1=trues_ind$value$alpha
            beta1=trues_ind$value$beta
            est_true_ind=c(est_true_ind,beta1[2])
            std_true_ind=c(std_true_ind,summary(trues_ind$value)$se.robust[2])
            w1=rep(1,dim(d1)[1])
            y1=d1$y
            X1=cbind(rep(1,dim(d1)[1]),d1$arm)
            id1=d1$cluster
            correction1=BCVexch(y1,X1,beta1,alpha1,phi1,id1,w1)
            robust_true_ind=c(robust_true_ind,sqrt(diag(correction1$robust))[2])
            MD_true_ind=c(MD_true_ind,sqrt(diag(correction1$varMD))[2])
            KC_true_ind=c(KC_true_ind,sqrt(diag(correction1$varKC))[2])
            FG_true_ind=c(FG_true_ind,sqrt(diag(correction1$varFG))[2])
          }else{
            est_true_ind=c(est_true_ind,NA)
            std_true_ind=c(std_true_ind,NA)
            robust_true_ind=c(robust_true_ind,NA)
            MD_true_ind=c(MD_true_ind,NA)
            KC_true_ind=c(KC_true_ind,NA)
            FG_true_ind=c(FG_true_ind,NA)
          }
        }else{
          est_true_ind=c(est_true_ind,NA)
          std_true_ind=c(std_true_ind,NA)
          robust_true_ind=c(robust_true_ind,NA)
          MD_true_ind=c(MD_true_ind,NA)
          KC_true_ind=c(KC_true_ind,NA)
          FG_true_ind=c(FG_true_ind,NA)
        }
        
        robust_true_ind;MD_true_ind;KC_true_ind;FG_true_ind
        est_true_ind;std_true_ind
        
        ## true exchangeable
        if(is.null(trues_ex$value)==0){
          if( trues_ex$value$converged==1){
            phi2=trues_ex$value$phi
            alpha2=trues_ex$value$alpha
            beta2=trues_ex$value$beta
            est_true_ex=c(est_true_ex,beta2[2])
            std_true_ex=c(std_true_ex,summary(trues_ex$value)$se.robust[2])
            w2=rep(1,dim(d1)[1])
            y2=d1$y
            X2=cbind(rep(1,dim(d1)[1]),d1$arm)
            id2=d1$cluster
            correction2=BCVexch(y2,X2,beta2,alpha2,phi2,id2,w2)
            robust_true_ex=c(robust_true_ex,sqrt(diag(correction2$robust))[2])
            MD_true_ex=c(MD_true_ex,sqrt(diag(correction2$varMD))[2])
            KC_true_ex=c(KC_true_ex,sqrt(diag(correction2$varKC))[2])
            FG_true_ex=c(FG_true_ex,sqrt(diag(correction2$varFG))[2])
          }else{
            est_true_ex=c(est_true_ex,NA)
            std_true_ex=c(std_true_ex,NA)
            robust_true_ex =c(robust_true_ex,NA)
            MD_true_ex=c(MD_true_ex,NA)
            KC_true_ex=c(KC_true_ex,NA)
            FG_true_ex=c(FG_true_ex,NA)
          }
        }else{
          est_true_ex=c(est_true_ex,NA)
          std_true_ex=c(std_true_ex,NA)
          robust_true_ex =c(robust_true_ex,NA)
          MD_true_ex=c(MD_true_ex,NA)
          KC_true_ex=c(KC_true_ex,NA)
          FG_true_ex=c(FG_true_ex,NA)
        }
        
        robust_true_ex;MD_true_ex;KC_true_ex;FG_true_ex
        est_true_ex;std_true_ex
        
        ### UCRA independent
        if(is.null(ucra_ind$value)==0){
          if( ucra_ind$value$converged==1){
            phi3=ucra_ind$value$phi
            alpha3=ucra_ind$value$alpha
            beta3=ucra_ind$value$beta
            est_ucra_ind=c(est_ucra_ind,beta3[2])
            std_ucra_ind=c(std_ucra_ind,summary(ucra_ind$value)$se.robust[2])
            w3=rep(1,dim(d2)[1])
            y3=d2$y
            X3=cbind(rep(1,dim(d2)[1]),d2$arm)
            id3=d2$cluster
            correction3=BCVexch(y3,X3,beta3,alpha3,phi3,id3,w3)
            robust_ucra_ind=c(robust_ucra_ind,sqrt(diag(correction3$robust))[2])
            MD_ucra_ind=c(MD_ucra_ind,sqrt(diag(correction3$varMD))[2])
            KC_ucra_ind=c(KC_ucra_ind,sqrt(diag(correction3$varKC))[2])
            FG_ucra_ind=c(FG_ucra_ind,sqrt(diag(correction3$varFG))[2])
          }else{
            est_ucra_ind=c(est_ucra_ind,NA)
            std_ucra_ind=c(std_ucra_ind,NA)
            robust_ucra_ind=c(robust_ucra_ind,NA)
            MD_ucra_ind=c(MD_ucra_ind,NA)
            KC_ucra_ind=c(KC_ucra_ind,NA)
            FG_ucra_ind=c(FG_ucra_ind,NA)
          }
        }else{
          est_ucra_ind=c(est_ucra_ind,NA)
          std_ucra_ind=c(std_ucra_ind,NA)
          robust_ucra_ind=c(robust_ucra_ind,NA)
          MD_ucra_ind=c(MD_ucra_ind,NA)
          KC_ucra_ind=c(KC_ucra_ind,NA)
          FG_ucra_ind=c(FG_ucra_ind,NA)
        }
        
        robust_ucra_ind;MD_ucra_ind;KC_ucra_ind;FG_ucra_ind
        est_ucra_ind;std_ucra_ind
        
        ## ucra exchangeable
        if(is.null(ucra_ex$value)==0 ){
          if(ucra_ex$value$converged==1){
            phi4=ucra_ex$value$phi
            alpha4=ucra_ex$value$alpha
            beta4=ucra_ex$value$beta
            est_ucra_ex=c(est_ucra_ex,beta4[2])
            std_ucra_ex=c(std_ucra_ex,summary(ucra_ex$value)$se.robust[2])
            w4=rep(1,dim(d2)[1])
            y4=d2$y
            X4=cbind(rep(1,dim(d2)[1]),d2$arm)
            id4=d2$cluster
            correction4=BCVexch(y4,X4,beta4,alpha4,phi4,id4,w4)
            robust_ucra_ex=c(robust_ucra_ex,sqrt(diag(correction4$robust))[2])
            MD_ucra_ex=c(MD_ucra_ex,sqrt(diag(correction4$varMD))[2])
            KC_ucra_ex=c(KC_ucra_ex,sqrt(diag(correction4$varKC))[2])
            FG_ucra_ex=c(FG_ucra_ex,sqrt(diag(correction4$varFG))[2])
            
          }else{
            est_ucra_ex=c(est_ucra_ex,NA)
            std_ucra_ex=c(std_ucra_ex,NA)
            robust_ucra_ex =c(robust_ucra_ex,NA)
            MD_ucra_ex=c(MD_ucra_ex,NA)
            KC_ucra_ex=c(KC_ucra_ex,NA)
            FG_ucra_ex=c(FG_ucra_ex,NA)
          }
        }else{
          est_ucra_ex=c(est_ucra_ex,NA)
          std_ucra_ex=c(std_ucra_ex,NA)
          robust_ucra_ex =c(robust_ucra_ex,NA)
          MD_ucra_ex=c(MD_ucra_ex,NA)
          KC_ucra_ex=c(KC_ucra_ex,NA)
          FG_ucra_ex=c(FG_ucra_ex,NA)
        }
        
        robust_ucra_ex;MD_ucra_ex;KC_ucra_ex;FG_ucra_ex
        est_ucra_ex;std_ucra_ex
        
        ## cra independent
        if(is.null(cra_ind$value)==0){
          if(cra_ind$value$converged==1){
            phi5=cra_ind$value$phi
            alpha5=cra_ind$value$alpha
            beta5=cra_ind$value$beta
            est_cra_ind=c(est_cra_ind,beta5[3])
            std_cra_ind=c(std_cra_ind,summary(cra_ind$value)$se.robust[2])
            w5=rep(1,dim(d2)[1])
            y5=d2$y
            X5=cbind(rep(1,dim(d2)[1]),d2$x,d2$arm)
            id5=d2$cluster
            correction5=BCVexch(y5,X5,beta5,alpha5,phi5,id5,w5)
            robust_cra_ind=c(robust_cra_ind,sqrt(diag(correction5$robust))[2])
            MD_cra_ind=c(MD_cra_ind,sqrt(diag(correction5$varMD))[2])
            KC_cra_ind=c(KC_cra_ind,sqrt(diag(correction5$varKC))[2])
            FG_cra_ind=c(FG_cra_ind,sqrt(diag(correction5$varFG))[2])
          }else{
            est_cra_ind=c(est_cra_ind,NA)
            std_cra_ind=c(std_cra_ind,NA)
            robust_cra_ind=c(robust_cra_ind,NA)
            MD_cra_ind=c(MD_cra_ind,NA)
            KC_cra_ind=c(KC_cra_ind,NA)
            FG_cra_ind=c(FG_cra_ind,NA)
          }
        }else{
          est_cra_ind=c(est_cra_ind,NA)
          std_cra_ind=c(std_cra_ind,NA)
          robust_cra_ind=c(robust_cra_ind,NA)
          MD_cra_ind=c(MD_cra_ind,NA)
          KC_cra_ind=c(KC_cra_ind,NA)
          FG_cra_ind=c(FG_cra_ind,NA)
        }
        robust_cra_ind;MD_cra_ind;KC_cra_ind;FG_cra_ind
        est_cra_ind;std_cra_ind
        
        ## cra exchangeable
        if(is.null(cra_ex$value)==0){
          if(cra_ex$value$converged==1){
            phi6=cra_ex$value$phi
            alpha6=cra_ex$value$alpha
            beta6=cra_ex$value$beta
            est_cra_ex=c(est_cra_ex,beta6[3])
            std_cra_ex=c(std_cra_ex,summary(cra_ex$value)$se.robust[2])
            w6=rep(1,dim(d1)[1])
            y6=d2$y
            X6=cbind(rep(1,dim(d2)[1]),d2$x,d2$arm)
            id6=d2$cluster
            correction6=BCVexch(y6,X6,beta6,alpha6,phi6,id6,w6)
            robust_cra_ex=c(robust_cra_ex,sqrt(diag(correction6$robust))[2])
            MD_cra_ex=c(MD_cra_ex,sqrt(diag(correction6$varMD))[2])
            KC_cra_ex=c(KC_cra_ex,sqrt(diag(correction6$varKC))[2])
            FG_cra_ex=c(FG_cra_ex,sqrt(diag(correction6$varFG))[2])
          }else{
            est_cra_ex=c(est_cra_ex,NA)
            std_cra_ex=c(std_cra_ex,NA)
            robust_cra_ex =c(robust_cra_ex,NA)
            MD_cra_ex=c(MD_cra_ex,NA)
            KC_cra_ex=c(KC_cra_ex,NA)
            FG_cra_ex=c(FG_cra_ex,NA)
          }
        }else{
          est_cra_ex=c(est_cra_ex,NA)
          std_cra_ex=c(std_cra_ex,NA)
          robust_cra_ex =c(robust_cra_ex,NA)
          MD_cra_ex=c(MD_cra_ex,NA)
          KC_cra_ex=c(KC_cra_ex,NA)
          FG_cra_ex=c(FG_cra_ex,NA)
        }
        robust_cra_ex;MD_cra_ex;KC_cra_ex;FG_cra_ex
        est_cra_ex;std_cra_ex
        
        ### IPW independent
        if(is.null(ipw_ind$value)==0){
          if(ipw_ind$value$converged==1){
            phi7=ipw_ind$value$phi
            alpha7=ipw_ind$value$alpha
            beta7=ipw_ind$value$beta
            est_ipw_ind=c(est_ipw_ind,beta7[2])
            std_ipw_ind=c( std_ipw_ind,summary(ipw_ind$value)$se.robust[2])
            w7=d4$weight
            y7=d4$y
            X7=cbind(rep(1,dim(d4)[1]),d4$arm)
            id7=d4$cluster
            correction7=BCVexch(y7,X7,beta7,alpha7,phi7,id7,w7)
            robust_ipw_ind=c(robust_ipw_ind,sqrt(diag(correction7$robust))[2])
            MD_ipw_ind=c(MD_ipw_ind,sqrt(diag(correction7$varMD))[2])
            KC_ipw_ind=c(KC_ipw_ind,sqrt(diag(correction7$varKC))[2])
            FG_ipw_ind=c(FG_ipw_ind,sqrt(diag(correction7$varFG))[2])
          }else{
            est_ipw_ind=c(est_ipw_ind,NA)
            std_ipw_ind=c(std_ipw_ind,NA)
            robust_ipw_ind=c(robust_ipw_ind,NA)
            MD_ipw_ind=c(MD_ipw_ind,NA)
            KC_ipw_ind=c(KC_ipw_ind,NA)
            FG_ipw_ind=c(FG_ipw_ind,NA)
          }
        }else{
          est_ipw_ind=c(est_ipw_ind,NA)
          std_ipw_ind=c(std_ipw_ind,NA)
          robust_ipw_ind=c(robust_ipw_ind,NA)
          MD_ipw_ind=c(MD_ipw_ind,NA)
          KC_ipw_ind=c(KC_ipw_ind,NA)
          FG_ipw_ind=c(FG_ipw_ind,NA)
        }
        
        robust_ipw_ind;MD_ipw_ind;KC_ipw_ind;FG_ipw_ind
        est_ipw_ind;std_ipw_ind
        
        ## ipw exchangeable
        if(is.null(ipw_ex$value)==0){
          if(ipw_ex$value$converged==1){
            phi8=ipw_ex$value$phi
            alpha8=ipw_ex$value$alpha
            beta8=ipw_ex$value$beta
            est_ipw_ex=c(est_ipw_ex,beta8[2])
            std_ipw_ex=c(std_ipw_ex,summary(ipw_ex$value)$se.robust[2])
            w8=d4$weight
            y8=d4$y
            X8=cbind(rep(1,dim(d4)[1]),d4$arm)
            id8=d4$cluster
            correction8=BCVexch(y8,X8,beta8,alpha8,phi8,id8,w8)
            robust_ipw_ex=c(robust_ipw_ex,sqrt(diag(correction8$robust))[2])
            MD_ipw_ex=c(MD_ipw_ex,sqrt(diag(correction8$varMD))[2])
            KC_ipw_ex=c(KC_ipw_ex,sqrt(diag(correction8$varKC))[2])
            FG_ipw_ex=c(FG_ipw_ex,sqrt(diag(correction8$varFG))[2])
          }else{
            est_ipw_ex=c(est_ipw_ex,NA)
            std_ipw_ex=c(std_ipw_ex,NA)
            robust_ipw_ex =c(robust_ipw_ex,NA)
            MD_ipw_ex=c(MD_ipw_ex,NA)
            KC_ipw_ex=c(KC_ipw_ex,NA)
            FG_ipw_ex=c(FG_ipw_ex,NA)
          }
        }else{
          est_ipw_ex=c(est_ipw_ex,NA)
          std_ipw_ex=c(std_ipw_ex,NA)
          robust_ipw_ex =c(robust_ipw_ex,NA)
          MD_ipw_ex=c(MD_ipw_ex,NA)
          KC_ipw_ex=c(KC_ipw_ex,NA)
          FG_ipw_ex=c(FG_ipw_ex,NA)
        }
        robust_ipw_ex;MD_ipw_ex;KC_ipw_ex;FG_ipw_ex
        est_ipw_ex;std_ipw_ex
        
        ### IPW_CLU independent
        if(is.null(ipw_clu_ind$value)==0){
          if(ipw_clu_ind$value$converged==1){
            phi9=ipw_clu_ind$value$phi
            alpha9=ipw_clu_ind$value$alpha
            beta9=ipw_clu_ind$value$beta
            est_ipw_clu_ind=c(est_ipw_clu_ind,beta9[2])
            std_ipw_clu_ind=c(std_ipw_clu_ind,summary(ipw_clu_ind$value)$se.robust[2])
            w9=d4$weight2
            y9=d4$y
            X9=cbind(rep(1,dim(d4)[1]),d4$arm)
            id9=d4$cluster
            correction9=BCVexch(y9,X9,beta9,alpha9,phi9,id9,w9)
            robust_ipw_clu_ind=c(robust_ipw_clu_ind,sqrt(diag(correction9$robust))[2])
            MD_ipw_clu_ind=c(MD_ipw_clu_ind,sqrt(diag(correction9$varMD))[2])
            KC_ipw_clu_ind=c(KC_ipw_clu_ind,sqrt(diag(correction9$varKC))[2])
            FG_ipw_clu_ind=c(FG_ipw_clu_ind,sqrt(diag(correction9$varFG))[2])
          }else{
            est_ipw_clu_ind=c(est_ipw_clu_ind,NA)
            std_ipw_clu_ind=c(std_ipw_clu_ind,NA)
            robust_ipw_clu_ind=c(robust_ipw_clu_ind,NA)
            MD_ipw_clu_ind=c(MD_ipw_clu_ind,NA)
            KC_ipw_clu_ind=c(KC_ipw_clu_ind,NA)
            FG_ipw_clu_ind=c(FG_ipw_clu_ind,NA)
          }
        }else{
          est_ipw_clu_ind=c(est_ipw_clu_ind,NA)
          std_ipw_clu_ind=c(std_ipw_clu_ind,NA)
          robust_ipw_clu_ind=c(robust_ipw_clu_ind,NA)
          MD_ipw_clu_ind=c(MD_ipw_clu_ind,NA)
          KC_ipw_clu_ind=c(KC_ipw_clu_ind,NA)
          FG_ipw_clu_ind=c(FG_ipw_clu_ind,NA)
        }
        
        robust_ipw_clu_ind;MD_ipw_clu_ind;KC_ipw_clu_ind;FG_ipw_clu_ind
        est_ipw_clu_ind;std_ipw_clu_ind
        
        ## ipw_clu exchangeable
        
        if(is.null(ipw_clu_ex$value)==0){
          if( ipw_clu_ex$value$converged==1){
            phi10=ipw_clu_ex$value$phi
            alpha10=ipw_clu_ex$value$alpha
            beta10=ipw_clu_ex$value$beta
            est_ipw_clu_ex=c(est_ipw_clu_ex,beta10[2])
            std_ipw_clu_ex=c(std_ipw_clu_ex,summary(ipw_clu_ex$value)$se.robust[2])
            w10=d4$weight2
            y10=d4$y
            X10=cbind(rep(1,dim(d4)[1]),d4$arm)
            id10=d4$cluster
            correction10=BCVexch(y10,X10,beta10,alpha10,phi10,id10,w10)
            robust_ipw_clu_ex=c(robust_ipw_clu_ex,sqrt(diag(correction10$robust))[2])
            MD_ipw_clu_ex=c(MD_ipw_clu_ex,sqrt(diag(correction10$varMD))[2])
            KC_ipw_clu_ex=c(KC_ipw_clu_ex,sqrt(diag(correction10$varKC))[2])
            FG_ipw_clu_ex=c(FG_ipw_clu_ex,sqrt(diag(correction10$varFG))[2])
          }else{
            est_ipw_clu_ex=c(est_ipw_clu_ex,NA)
            std_ipw_clu_ex=c(std_ipw_clu_ex,NA)
            robust_ipw_clu_ex =c(robust_ipw_clu_ex,NA)
            MD_ipw_clu_ex=c(MD_ipw_clu_ex,NA)
            KC_ipw_clu_ex=c(KC_ipw_clu_ex,NA)
            FG_ipw_clu_ex=c(FG_ipw_clu_ex,NA)
          }
        }else{
          est_ipw_clu_ex=c(est_ipw_clu_ex,NA)
          std_ipw_clu_ex=c(std_ipw_clu_ex,NA)
          robust_ipw_clu_ex =c(robust_ipw_clu_ex,NA)
          MD_ipw_clu_ex=c(MD_ipw_clu_ex,NA)
          KC_ipw_clu_ex=c(KC_ipw_clu_ex,NA)
          FG_ipw_clu_ex=c(FG_ipw_clu_ex,NA)
        }
        robust_ipw_clu_ex;MD_ipw_clu_ex;KC_ipw_clu_ex;FG_ipw_clu_ex
        est_ipw_clu_ex;std_ipw_clu_ex
      }
      
      est_ind=data.frame(est_true_ind=est_true_ind,est_ucra_ind=est_ucra_ind,
                         est_cra_ind=est_cra_ind,est_ipw_ind=est_ipw_ind,
                         est_ipw_clu_ind=est_ipw_clu_ind)
      est_ex=data.frame(est_true_ex=est_true_ex,est_ucra_ex=est_ucra_ex,
                        est_cra_ex=est_cra_ex,est_ipw_ex=est_ipw_ex,
                        est_ipw_clu_ex=est_ipw_clu_ex)
      
      std_ind=data.frame(std_true_ind=std_true_ind,std_ucra_ind=std_ucra_ind,
                         std_cra_ind=std_cra_ind,std_ipw_ind=std_ipw_ind,
                         std_ipw_clu_ind=std_ipw_clu_ind)
      std_ex=data.frame(std_true_ex=std_true_ex,std_ucra_ex=std_ucra_ex,
                        std_cra_ex=std_cra_ex,std_ipw_ex=std_ipw_ex,
                        std_ipw_clu_ex=std_ipw_clu_ex)
      
      robust_ind=data.frame(robust_true_ind=robust_true_ind,robust_ucra_ind=robust_ucra_ind,
                            robust_cra_ind=robust_cra_ind,robust_ipw_ind=robust_ipw_ind,
                            robust_ipw_clu_ind=robust_ipw_clu_ind)
      robust_ex=data.frame(robust_true_ex=robust_true_ex,robust_ucra_ex=robust_ucra_ex,
                           robust_cra_ex=robust_cra_ex,robust_ipw_ex=robust_ipw_ex,
                           robust_ipw_clu_ex=robust_ipw_clu_ex)
      
      MD_ind=data.frame(MD_true_ind=MD_true_ind,MD_ucra_ind=MD_ucra_ind,
                        MD_cra_ind=MD_cra_ind,MD_ipw_ind=MD_ipw_ind,
                        MD_ipw_clu_ind=MD_ipw_clu_ind)
      MD_ex=data.frame(MD_true_ex=MD_true_ex,MD_ucra_ex=MD_ucra_ex,
                       MD_cra_ex=MD_cra_ex,MD_ipw_ex=MD_ipw_ex,
                       MD_ipw_clu_ex=MD_ipw_clu_ex)
      
      
      KC_ind=data.frame(KC_true_ind=KC_true_ind,KC_ucra_ind=KC_ucra_ind,
                        KC_cra_ind=KC_cra_ind,KC_ipw_ind=KC_ipw_ind,
                        KC_ipw_clu_ind=KC_ipw_clu_ind)
      KC_ex=data.frame(KC_true_ex=KC_true_ex,KC_ucra_ex=KC_ucra_ex,
                       KC_cra_ex=KC_cra_ex,KC_ipw_ex=KC_ipw_ex,
                       KC_ipw_clu_ex=KC_ipw_clu_ex)
      
      FG_ind=data.frame(FG_true_ind=FG_true_ind,FG_ucra_ind=FG_ucra_ind,
                        FG_cra_ind=FG_cra_ind,FG_ipw_ind=FG_ipw_ind,
                        FG_ipw_clu_ind=FG_ipw_clu_ind)
      FG_ex=data.frame(FG_true_ex=FG_true_ex,FG_ucra_ex=FG_ucra_ex,
                       FG_cra_ex=FG_cra_ex,FG_ipw_ex=FG_ipw_ex,
                       FG_ipw_clu_ex=FG_ipw_clu_ex)
      
      result=list(est_ind=est_ind,est_ex=est_ex,std_ind=std_ind,std_ex=std_ex,
                  robust_ind=robust_ind,robust_ex=robust_ex,
                  MD_ind=MD_ind,KC_ind=KC_ind,FG_ind=FG_ind,
                  MD_ex=MD_ex,KC_ex=KC_ex,FG_ex=FG_ex)
      names=paste('geem',k,icc,iccm,'.RData',sep='')
      save(result,file=names)
    }
  }
  
}

### 5. Combine results and make table

## draw ipw cra tables
setwd('')

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
      biasind=100*(estind[2:5]-estind[1])/estind[1]  ## bias
      biasex=100*(estex[2:5]-estex[1])/estex[1]   ## bias
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
      setwd('')
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
      cov_mmi_ind=sum(((est_ind-1.96*std_ind)<trueind) &  ((est_ind+1.96*std_ind)>trueind))/1000
      cov_mmi_ex=sum(((est_ex-1.96*std_ex)<trueex) &  ((est_ex+1.96*std_ex)>trueex))/1000
      
      cov_mmi_ro_ind=sum(((est_ind-1.96*robust_ind)<trueind) &  ((est_ind+1.96*robust_ind)>trueind))/1000
      cov_mmi_ro_ex=sum(((est_ex-1.96*robust_ex)<trueex) &  ((est_ex+1.96*robust_ex)>trueex))/1000
      
      cov_mmi_MD_ind=sum(((est_ind-1.96*md_ind)<trueind) &  ((est_ind+1.96*md_ind)>trueind))/1000
      cov_mmi_MD_ex=sum(((est_ex-1.96*md_ex)<trueex) &  ((est_ex+1.96*md_ex)>trueex))/1000
      
      cov_mmi_KC_ind=sum(((est_ind-1.96*kc_ind)<trueind) &  ((est_ind+1.96*kc_ind)>trueind))/1000
      cov_mmi_KC_ex=sum(((est_ex-1.96*kc_ex)<trueex) &  ((est_ex+1.96*kc_ex)>trueex))/1000
      
      cov_mmi_FG_ind=sum(((est_ind-1.96*fg_ind)<trueind) &  ((est_ind+1.96*fg_ind)>trueind))/1000
      cov_mmi_FG_ex=sum(((est_ex-1.96*fg_ex)<trueex) &  ((est_ex+1.96*fg_ex)>trueex))/1000
      
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

## FINDAL TABLES 
table_ind=cbind(BIAS_ind,Cov_ind[,4:8],STD_ind[,4:8],MCSD_ind[,4:8],uncov_ind[,5:9])
table_ex=cbind(BIAS_ex,Cov_ex[,4:8],STD_ex[,4:8],MCSD_ex[,4:8],uncov_ex[,5:9])

