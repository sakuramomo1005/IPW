
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


### simulation

for(icc in c(0.01,0.05,0.1)){
true_est_ind=c();true_std_ind=c();true_warn_ind=c()
true_est_ex=c();true_std_ex=c();true_warn_ex=c()
ucra_est_ind=c();ucra_std_ind=c();ucra_warn_ind=c()
ucra_est_ex=c();ucra_std_ex=c();ucra_warn_ex=c()
cra_est_ind=c();cra_std_ind=c();cra_warn_ind=c()
cra_est_ex=c();cra_std_ex=c();cra_warn_ex=c()
ipw_est_ind=c();ipw_std_ind=c();ipw_warn_ind=c()
ipw_est_ex=c();ipw_std_ex=c();ipw_warn_ex=c()
ipw_clu_est_ind=c();ipw_clu_std_ind=c();ipw_clu_warn_ind=c()
ipw_clu_est_ex=c();ipw_clu_std_ex=c();ipw_clu_warn_ex=c()

mis=c()

for(times in 1:1000){
  
  print(times)
  vd=missing_icc(icc)
  
  d1=datagen(vx=2,vd=vd,psi=-1.6,seed=times)
  mis=c(mis,sum(d1$missing))
  
  for(i in 1:100){
    t=d1[d1$cluster==i,]
    if(sum(t$missing)==50) break
  }
  
  d2=na.omit(d1)
  d3=d1
  d3$y=ifelse(d3$missing==1,NA,d3$y)
  
  w1=glm(missing ~ x+arm, data = d3,
         family = binomial(link='logit'))
  w2=glmer(missing ~ x+arm+(1|cluster) , data = d3,
           family = binomial(link='logit'))
  
  weight1=1/predict(w1,type="response")
  weight2=1/predict(w2,type="response")
  
  d3$weight=weight1
  d3$weight2=weight2
  
  ### True effect
  trues_ind=myTryCatch(geeDREstimation(formula=y~x+arm,
                                       id="cluster" , data = d1,nameY='y',
                                       nameTRT='arm',
                                       family =  binomial("logit"),
                                       corstr = "independence"))
  trues_ex=myTryCatch(geeDREstimation(formula=y~x+arm,
                                      id="cluster" , data = d1,nameY='y',
                                      nameTRT='arm',
                                      family =  binomial("logit"),
                                      corstr = "exchangeable"))
  
  ### Unadjusted CRA
  ucra_ind=myTryCatch(geeDREstimation(formula=y~arm,
                                      id="cluster" , data = d2,
                                      nameMISS='missing',nameY='y',
                                      nameTRT='arm',
                                      family =  binomial("logit"),
                                      corstr = "independence"))
  ucra_ex=myTryCatch(geeDREstimation(formula=y~arm,
                                     id="cluster" , data = d2,
                                     nameMISS='missing',nameY='y',
                                     nameTRT='arm',
                                     family =  binomial("logit"),
                                     corstr = "exchangeable"))
  
  ### Adjusted CRA
  cra_ind=myTryCatch(geeDREstimation(formula=y~x+arm,
                                     id="cluster" , data = d2,
                                     nameMISS='missing',nameY='y',
                                     nameTRT='arm',
                                     family =  binomial("logit"),
                                     corstr = "independence"))
  cra_ex=myTryCatch(geeDREstimation(formula=y~x+arm,
                                    id="cluster" , data = d2,
                                    nameMISS='missing',nameY='y',
                                    nameTRT='arm',
                                    family =  binomial("logit"),
                                    corstr = "exchangeable"))
  
  ### IPW
  ipw_ind=myTryCatch(geeDREstimation(formula=y~x+arm,
                                     id="cluster" , data = d3,
                                     nameMISS='missing',nameY='y',
                                     nameTRT='arm',
                                     weights = d3$weight,
                                     family =  binomial("logit"),
                                     corstr = "independence"))
  
  ipw_ex=myTryCatch(geeDREstimation(formula=y~x+arm,
                                    id="cluster" , data = d3,
                                    nameMISS='missing',nameY='y',
                                    nameTRT='arm',
                                    weights = d3$weight,
                                    family =  binomial("logit"),
                                    corstr = "exchangeable"))
  
  ### IPW-CLU
  ipw_clu_ind=myTryCatch(geeDREstimation(formula=y~x+arm,
                                         id="cluster" , data = d3,
                                         nameMISS='missing',nameY='y',
                                         nameTRT='arm',
                                         weights = d3$weight2,
                                         family =  binomial("logit"),
                                         corstr = "independence"))
  
  ipw_clu_ex=myTryCatch(geeDREstimation(formula=y~x+arm,
                                        id="cluster" , data = d3,
                                        nameMISS='missing',nameY='y',
                                        nameTRT='arm',
                                        weights = d3$weight2,
                                        family =  binomial("logit"),
                                        corstr = "exchangeable"))
  
  #### save results:
  
  
  ## true
  if(is.null(trues_ind$value)==0){
    t1=summary(trues_ind$value)$beta[3]
    t2=summary(trues_ind$value)$se.robust[3]
    true_est_ind=c(true_est_ind,t1)
    true_std_ind=c(true_std_ind,t2)
    if(trues_ind$value$converged==0){true_warn_ind=c(true_warn_ind,times)}
    if(trues_ind$value$converged==1){true_warn_ind=c(true_warn_ind,0)}
  }
  if(is.null(trues_ex$value)==0){
    t1=summary(trues_ex$value)$beta[3]
    t2=summary(trues_ex$value)$se.robust[3]
    true_est_ex=c(true_est_ex,t1)
    true_std_ex=c(true_std_ex,t2)
    if(trues_ex$value$converged==0){true_warn_ex=c(true_warn_ex,times)}
    if(trues_ex$value$converged==1){true_warn_ex=c(true_warn_ex,0)}
  }
  
  if(is.null(trues_ind$value)==1){
    true_est_ind=c(true_est_ind,NA)
    true_std_ind=c(true_std_ind,NA)
    true_warn_ind=c(true_warn_ind,times)
  }
  if(is.null(trues_ex$value)==1){
    true_est_ex=c(true_est_ex,NA)
    true_std_ex=c(true_std_ex,NA)
    true_warn_ex=c(true_warn_ex,times)
  }
  
  ### UCRA
  if(is.null(ucra_ind$value)==0){
    t1=summary(ucra_ind$value)$beta[2]
    t2=summary(ucra_ind$value)$se.robust[2]
    ucra_est_ind=c(ucra_est_ind,t1)
    ucra_std_ind=c(ucra_std_ind,t2)
    if(ucra_ind$value$converged==0){ucra_warn_ind=c(ucra_warn_ind,times)}
    if(ucra_ind$value$converged==1){ucra_warn_ind=c(ucra_warn_ind,0)}
  }
  if(is.null(ucra_ex$value)==0){
    t1=summary(ucra_ex$value)$beta[2]
    t2=summary(ucra_ex$value)$se.robust[2]
    ucra_est_ex=c(ucra_est_ex,t1)
    ucra_std_ex=c(ucra_std_ex,t2)
    if(ucra_ex$value$converged==0){ucra_warn_ex=c(ucra_warn_ex,times)}
    if(ucra_ex$value$converged==1){ucra_warn_ex=c(ucra_warn_ex,0)}
  }
  
  if(is.null(ucra_ind$value)==1){
    ucra_est_ind=c(ucra_est_ind,NA)
    ucra_std_ind=c(ucra_std_ind,NA)
    ucra_warn_ind=c(ucra_warn_ind,times)
  }
  if(is.null(ucra_ex$value)==1){
    ucra_est_ex=c(ucra_est_ex,NA)
    ucra_std_ex=c(ucra_std_ex,NA)
    ucra_warn_ex=c(ucra_warn_ex,times)
  }
  
  
  ### CRA
  if(is.null(cra_ind$value)==0){
    t1=summary(cra_ind$value)$beta[3]
    t2=summary(cra_ind$value)$se.robust[3]
    cra_est_ind=c(cra_est_ind,t1)
    cra_std_ind=c(cra_std_ind,t2)
    if(cra_ind$value$converged==0){cra_warn_ind=c(cra_warn_ind,times)}
    if(cra_ind$value$converged==1){cra_warn_ind=c(cra_warn_ind,0)}
  }
  if(is.null(cra_ex$value)==0){
    t1=summary(cra_ex$value)$beta[3]
    t2=summary(cra_ex$value)$se.robust[3]
    cra_est_ex=c(cra_est_ex,t1)
    cra_std_ex=c(cra_std_ex,t2)
    if(cra_ex$value$converged==0){cra_warn_ex=c(cra_warn_ex,times)}
    if(cra_ex$value$converged==1){cra_warn_ex=c(cra_warn_ex,0)}
  }
  
  if(is.null(cra_ind$value)==1){
    cra_est_ind=c(cra_est_ind,NA)
    cra_std_ind=c(cra_std_ind,NA)
    cra_warn_ind=c(cra_warn_ind,times)
  }
  if(is.null(cra_ex$value)==1){
    cra_est_ex=c(cra_est_ex,NA)
    cra_std_ex=c(cra_std_ex,NA)
    cra_warn_ex=c(cra_warn_ex,times)
  }
  
  
  ### IPW
  if(is.null(ipw_ind$value)==0){
    t1=summary(ipw_ind$value)$beta[3]
    t2=summary(ipw_ind$value)$se.robust[3]
    ipw_est_ind=c(ipw_est_ind,t1)
    ipw_std_ind=c(ipw_std_ind,t2)
    if(ipw_ind$value$converged==0){ipw_warn_ind=c(ipw_warn_ind,times)}
    if(ipw_ind$value$converged==1){ipw_warn_ind=c(ipw_warn_ind,0)}
  }
  if(is.null(ipw_ex$value)==0){
    t1=summary(ipw_ex$value)$beta[3]
    t2=summary(ipw_ex$value)$se.robust[3]
    ipw_est_ex=c(ipw_est_ex,t1)
    ipw_std_ex=c(ipw_std_ex,t2)
    if(ipw_ex$value$converged==0){ipw_warn_ex=c(ipw_warn_ex,times)}
    if(ipw_ex$value$converged==1){ipw_warn_ex=c(ipw_warn_ex,0)}
  }
  
  if(is.null(ipw_ind$value)==1){
    ipw_est_ind=c(ipw_est_ind,NA)
    ipw_std_ind=c(ipw_std_ind,NA)
    ipw_warn_ind=c(ipw_warn_ind,times)
  }
  if(is.null(ipw_ex$value)==1){
    ipw_est_ex=c(ipw_est_ex,NA)
    ipw_std_ex=c(ipw_std_ex,NA)
    ipw_warn_ex=c(ipw_warn_ex,times)
  }
  
  ### IPW cluster
  if(is.null(ipw_clu_ind$value)==0){
    t1=summary(ipw_clu_ind$value)$beta[3]
    t2=summary(ipw_clu_ind$value)$se.robust[3]
    ipw_clu_est_ind=c(ipw_clu_est_ind,t1)
    ipw_clu_std_ind=c(ipw_clu_std_ind,t2)
    if(ipw_clu_ind$value$converged==0){ipw_clu_warn_ind=c(ipw_clu_warn_ind,times)}
    if(ipw_clu_ind$value$converged==1){ipw_clu_warn_ind=c(ipw_clu_warn_ind,0)}
  }
  if(is.null(ipw_clu_ex$value)==0){
    t1=summary(ipw_clu_ex$value)$beta[3]
    t2=summary(ipw_clu_ex$value)$se.robust[3]
    ipw_clu_est_ex=c(ipw_clu_est_ex,t1)
    ipw_clu_std_ex=c(ipw_clu_std_ex,t2)
    if(ipw_clu_ex$value$converged==0){ipw_clu_warn_ex=c(ipw_clu_warn_ex,times)}
    if(ipw_clu_ex$value$converged==1){ipw_clu_warn_ex=c(ipw_clu_warn_ex,0)}
  }
  if(is.null(ipw_clu_ind$value)==1){
    ipw_clu_est_ind=c(ipw_clu_est_ind,NA)
    ipw_clu_std_ind=c(ipw_clu_std_ind,NA)
    ipw_clu_warn_ind=c(ipw_clu_warn_ind,times)
  }
  if(is.null(ipw_clu_ex$value)==1){
    ipw_clu_est_ex=c(ipw_clu_est_ex,NA)
    ipw_clu_std_ex=c(ipw_clu_std_ex,NA)
    ipw_clu_warn_ex=c(ipw_clu_warn_ex,times)
  }
  
}


est_ind=data.frame(true_est_ind=true_est_ind,ucra_est_ind=ucra_est_ind,cra_est_ind=cra_est_ind,
                   ipw_est_ind=ipw_est_ind,ipw_clu_est_ind=ipw_clu_est_ind)
est_ex=data.frame(true_est_ex=true_est_ex,ucra_est_ex=ucra_est_ex,cra_est_ex=cra_est_ex,
                  ipw_est_ex=ipw_est_ex,ipw_clu_est_ex=ipw_clu_est_ex)

std_ind=data.frame(true_std_ind=true_std_ind,ucra_std_ind=ucra_std_ind,cra_std_ind=cra_std_ind,
                   ipw_std_ind=ipw_std_ind,ipw_clu_std_ind=ipw_clu_std_ind)
std_ex=data.frame(true_std_ex=true_std_ex,ucra_std_ex=ucra_std_ex,cra_std_ex=cra_std_ex,
                  ipw_std_ex=ipw_std_ex,ipw_clu_std_ex=ipw_clu_std_ex)

warn_ind=data.frame(true_warn_ind=true_warn_ind,ucra_warn_ind=ucra_warn_ind,cra_warn_ind=cra_warn_ind,
                    ipw_warn_ind=ipw_warn_ind,ipw_clu_warn_ind=ipw_clu_warn_ind)
warn_ex=data.frame(true_warn_ex=true_warn_ex,ucra_warn_ex=ucra_warn_ex,cra_warn_ex=cra_warn_ex,
                   ipw_warn_ex=ipw_warn_ex,ipw_clu_warn_ex=ipw_clu_warn_ex)

result=list(est_ind=est_ind,est_ex=est_ex,std_ind=std_ind,
            std_ex=std_ex,warn_ind,warn_ind,warn_ex=warn_ex)
names=paste('res',icc,'.RData',sep='')
save(result,file=names)
}