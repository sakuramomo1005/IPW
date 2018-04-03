library(jomo)
library(lme4)
library(CRTgeeDR)


## function to catch errors and warns
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

# expit function
expit=function(x){y=exp(x)/(1+exp(x));return(y)}

## from missingness ICC calculate variance
missing_icc=function(rho){
  pi=3.1415926
  sigma=(pi^2/3)/(1/rho-1)
  return(sigma)
}

## data generation
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

## pool function
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


Nimp=15

for(icc in c(0.01,0.05,0.1)){
  est_ind=c();est_ex=c()
  std_ind=c();std_ex=c()
  warn_ind=c();warn_ex=c()
  
for(times in 1:1000){
  print('times')
  print(times)
  vd=missing_icc(icc)
  d1=datagen(vx=2,vd=vd,psi=-1.6,seed=times)
  d3=d1
  d3$y=ifelse(d3$missing==1,NA,d3$y)
  data.miss=d3
  
  y.cat= data.frame(outcome=data.miss$y)  # data frame for response variables with missing values
  y.numcat=c(2)                                 # number of levels in outcome variable
  clus=data.frame(clus=data.miss$cluster)          # data frame for clusters
  
  nobs=dim(data.miss)[1]
  x= data.frame(intercept=rep(1,nobs),covariate=data.miss$x,group=data.miss$arm) 
  imp = jomo1rancat(Y.cat=y.cat, Y.numcat=y.numcat, X=x,
                    clus=clus,nburn=100, nbetween=25, nimp=Nimp,output=0)
  
  
  mmi_est_ind=c();mmi_std_ind=c();mmi_warn_ind=c()
  mmi_est_ex=c();mmi_std_ex=c();mmi_warn_ex=c()
  for(i in 1:Nimp){
    temp=imp[imp$Imputation==i,]
    rownames(temp)=NULL
    temp$outcome=as.numeric(temp$outcome)-1
    
    mmi_ind=myTryCatch(geeDREstimation(formula=outcome~covariate+group,
                                       id="clus" , data = temp, nameY='outcome',
                                       nameTRT='group',
                                       family =  binomial("logit"),
                                       corstr = "independence"))
    
    mmi_ex=myTryCatch(geeDREstimation(formula=outcome~covariate+group,
                                      id="clus" , data = temp, nameY='outcome',
                                      nameTRT='group',
                                      family =  binomial("logit"),
                                      corstr = "exchangeable"))
    ## mmi
    if(is.null(mmi_ind$value)==0){
      t1=summary(mmi_ind$value)$beta[3]
      t2=summary(mmi_ind$value)$se.robust[3]
      mmi_est_ind=c(mmi_est_ind,t1)
      mmi_std_ind=c(mmi_std_ind,t2)
      if(mmi_ind$value$converged==0){mmi_warn_ind=c(mmi_warn_ind,times)}
      if(mmi_ind$value$converged==1){mmi_warn_ind=c(mmi_warn_ind,0)}
    }
    if(is.null(mmi_ex$value)==0){
      t1=summary(mmi_ex$value)$beta[3]
      t2=summary(mmi_ex$value)$se.robust[3]
      mmi_est_ex=c(mmi_est_ex,t1)
      mmi_std_ex=c(mmi_std_ex,t2)
      if(mmi_ex$value$converged==0){mmi_warn_ex=c(mmi_warn_ex,times)}
      if(mmi_ex$value$converged==1){mmi_warn_ex=c(mmi_warn_ex,0)}
    }
    
    if(is.null(mmi_ind$value)==1){
      mmi_est_ind=c(mmi_est_ind,NA)
      mmi_std_ind=c(mmi_std_ind,NA)
      mmi_warn_ind=c(mmi_warn_ind,times)
    }
    if(is.null(mmi_ex$value)==1){
      mmi_est_ex=c(mmi_est_ex,NA)
      mmi_std_ex=c(mmi_std_ex,NA)
      mmi_warn_ex=c(mmi_warn_ex,times)
    }
  }
  
  temp1=mypool(mmi_est_ind,mmi_std_ind,num=Nimp)
  temp2=mypool(mmi_est_ex,mmi_std_ex,num=Nimp)
  
  est_ind=c(est_ind,temp1$mean)
  std_ind=c(std_ind,temp1$std)
  est_ex=c(est_ex,temp2$mean)
  std_ex=c(std_ex,temp2$std)
  warn_ind=data.frame(mmi_warn_ind=mmi_warn_ind,n=i)
  warn_ex=data.frame(mmi_warn_ex=mmi_warn_ex,n=i)

}

result=list(est_ind=est_ind,est_ex=est_ex,std_ind=std_ind,
            std_ex=std_ex,warn_ind,warn_ind,warn_ex=warn_ex)
names=paste('res_mmi',icc,'.RData',sep='')
save(result,file=names)
}

