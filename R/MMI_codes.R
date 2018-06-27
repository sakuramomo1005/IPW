#####################################################################################
#                                    Roadmap
# This is the R code for:
#  1. Generate simulation datasets (datasets with missing outcomes)
#  2. Analyze the datasets with Multilevel Multiple Imputation (MMI)

# Details:
# The simulation steps are:
#  1. Generate a full datasets without missingness
#  2. Generate missing values in outcome based on the function we mention in the paper
#  3. Analyze the datasets with MMI
#  4. Save the results (mean values, sd, mcsd, coverage rate, non-convergence times) in
#     a RData file
#  5. Repeat step 1-4 for 1000 times
#  6. Load in the saved RData files to make a summary table

# MMI
#  1. Impute the dataset for 15 times
#  2. Analyze each imputated dataset with GEE function
#  3. Pool all the results
#####################################################################################


# load in the library. Use the MMI.  
library(jomo)

# Function to catch errors and warns
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

# Expit function
expit=function(x){y=exp(x)/(1+exp(x));return(y)}

# From missingness ICC calculate variance
missing_icc=function(rho){
  pi=3.1415926
  sigma=(pi^2/3)/(1/rho-1)
  return(sigma)
}

# Data generation
datagen=function(vx=1,vd,psi=-1.7,seed=123){
  set.seed(seed)
  cluster=rep(1:100,each=50)
  cluster=factor(cluster)
  arm=rep(0:1,each=2500)
  x=rnorm(5000,0,sqrt(vx))
  delta=rep(rnorm(100,0,sqrt(vd)),each=50)
  pi=expit(1+1.36*arm+x+delta)
  r=expit(psi+arm+x)
  y=c();missing=c()
  for(i in 1:5000){
    y=c(y,rbinom(1,1,pi[i]))
    missing=c(missing,rbinom(1,1,r[i]))
  }
  res=data.frame(x=x,y=y,missing=missing,pi=pi,arm=arm,cluster=cluster)
  return(res)
}

# Pool function
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



s=Sys.time()

T=1000

for(k in c(10, 25, 50)){
  for(icc in c(0.05,0.1,0.2)){
    for(iccm in c(0,0.1,0.3,0.5,0.7,0.9)){
      
      # set empty variables to store the results
      est_ind=c();est_ex=c()
      std_ind=c();std_ex=c()
      warn_ind=c();warn_ex=c()
      
      for(times in 1:T){
        print(paste('icc',icc,'iccm',iccm,'times',times))
        set.seed(times)
        
        ## Tune the intercept
        # since the parameters change may bring different missing percentages,
        # we tune the intercept to keep the missing percentage around 30%
        
        mis=c()
        for(intercept in seq(-5,5,0.1)){
          temp=dategen(k,M,varx=varx,icc=icc,iccm=iccm,intercept=intercept) 
          mis=c(mis,missing_per(temp))
        }
        intercept=seq(-5,5,0.1)[which.min(abs(mis-0.3))]
        
        
        # generate a dataset based on the tuned intercept
        # d1: the full dataset
        # d3: the dataset with misisng outcomes
        
        d1=dategen(k,M,varx=varx,icc=icc,iccm=iccm,intercept=intercept)
        d3=d1
        d3$y=ifelse(d3$r==1,NA,d3$y)
        d3$missing=d3$r
        
        data.miss=d3
        
        y.cat= data.frame(outcome=data.miss$y)  # data frame for response variables with missing values
        y.numcat=c(2)                                 # number of levels in outcome variable
        clus=data.frame(clus=data.miss$cluster)          # data frame for clusters
        
        nobs=dim(data.miss)[1]
        x= data.frame(intercept=rep(1,nobs),covariate=data.miss$x,group=data.miss$arm)
        
        # run to generate Nimp full datasets
        # imp is full datasets 
        imp = jomo1rancat(Y.cat=y.cat, Y.numcat=y.numcat, X=x,
                          clus=clus,nburn=100, nbetween=25, nimp=Nimp,output=0)
        mmi_est_ind=c();mmi_std_ind=c();mmi_warn_ind=c()
        mmi_est_ex=c();mmi_std_ex=c();mmi_warn_ex=c()
        
        # Analyze each of the full dataset with GEE
        for(i in 1:Nimp){
          temp=imp[imp$Imputation==i,]
          rownames(temp)=NULL
          temp$outcome=as.numeric(temp$outcome)-1
          
          mmi_ind=myTryCatch(geeglm(formula=outcome~group,
                                    id=clus , data = temp,
                                    family =  binomial("logit"),
                                    corstr = "independence"))
          
          mmi_ex=myTryCatch(geeglm(formula=outcome~group,
                                   id=clus , data = temp,
                                   family =  binomial("logit"),
                                   corstr = "exchangeable"))
          ## Save the results
          if(is.null(mmi_ind$value)==0){
            t1=summary(mmi_ind$value)$coefficients['group','Estimate']
            t2=summary(mmi_ind$value)$coefficients['group','Std.err']
            mmi_est_ind=c(mmi_est_ind,t1)
            mmi_std_ind=c(mmi_std_ind,t2)
            if(summary(mmi_ind$value)$error==1){mmi_warn_ind=c(mmi_warn_ind,times)}
            if(summary(mmi_ind$value)$error==0){mmi_warn_ind=c(mmi_warn_ind,0)}
          }
          if(is.null(mmi_ind$value)==1){
            mmi_est_ind=c(mmi_est_ind,NA)
            mmi_std_ind=c(mmi_std_ind,NA)
            mmi_warn_ind=c(mmi_warn_ind,times)
          }
          
          if(is.null(mmi_ex$value)==0){
            t1=summary(mmi_ex$value)$coefficients['group','Estimate']
            t2=summary(mmi_ex$value)$coefficients['group','Std.err']
            mmi_est_ex=c(mmi_est_ex,t1)
            mmi_std_ex=c(mmi_std_ex,t2)
            if(summary(mmi_ex$value)$error==1){mmi_warn_ex=c(mmi_warn_ex,times)}
            if(summary(mmi_ex$value)$error==0){mmi_warn_ex=c(mmi_warn_ex,0)}
          }
          if(is.null(mmi_ex$value)==1){
            mmi_est_ex=c(mmi_est_ex,NA)
            mmi_std_ex=c(mmi_std_ex,NA)
            mmi_warn_ex=c(mmi_warn_ex,times)
          }
          
        }
        
        # Pool the results together
        temp1=mypool(mmi_est_ind,mmi_std_ind,num=Nimp)
        temp2=mypool(mmi_est_ex,mmi_std_ex,num=Nimp)
        
        est_ind=c(est_ind,temp1$mean)
        std_ind=c(std_ind,temp1$std)
        est_ex=c(est_ex,temp2$mean)
        std_ex=c(std_ex,temp2$std)
        warn_ind=c(warn_ind,sum(mmi_warn_ind))
        warn_ex=c(warn_ex,sum(mmi_warn_ex))
      }
      
      result=list(est_ind=est_ind,est_ex=est_ex,std_ind=std_ind,
                  std_ex=std_ex,warn_ind,warn_ind,warn_ex=warn_ex)
      print(result)
      names=paste('wrapk_mmi',k,icc,iccm,'.RData',sep='')
      save(result,file=names)
    }
  }
}

e=Sys.time()
