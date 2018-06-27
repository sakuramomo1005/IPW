
#####################################################################################
#                                    Roadmap
#This is the R code for:
#  1. Generate simulation datasets (datasets with missing outcomes)
#  2. Analyze the datasets with unadjusted complete records analysis, adjusted complete
#     records analysis, inverse probability weighting without cluster effets, and
#     inverse probaility weighting with cluster effects

# Details:
# The simulation steps are:
#  1. Generate a full datasets without missingness
#  2. Generate missing values in outcome based on the function we mention in the paper
#  3. Analyze the datasets with above methods
#  4. Save the results (mean values, sd, mcsd, coverage rate, non-convergence times) in
#     a RData file
#  5. Repeat step 1-4 for 1000 times
#  6. Load in the saved RData files to make a summary table
#####################################################################################


# load in the needed libraries
library(lme4)
library(geepack)


#####################################################################################
#                                    Functions
#####################################################################################

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

# Function to calculate the variance based on ICC
missing_icc=function(icc){
  pi=3.142
  a=pi^2/(3*(1/icc-1))
  return(a)
}

# Function to calculate the missingness percentage in dataset
missing_per=function(data){
  res=sum(data$r)/dim(data)[1]
  return(res)
}

# Expit function
expit=function(x){y=exp(x)/(1+exp(x));return(y)}

# Data generation function
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


#####################################################################################
#                               Missing data analysis
#####################################################################################

setwd()  # Choose a file path to save the RData files

T=1000

s1=Sys.time()

for(k in c(10, 25, 50)){  # k, the number of clusters in each arm
  for(icc in c(0.05,0.1,0.2)){  # icc, the dataset's icc
    for(iccm in c(0, 0.1, 0.3, 0.5, 0.7, 0.9)){ # the missingness icc
      
      # set empty variables to save values
      
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
      miss=c()
        
      for(times in 1:T){ # T=1000, the repeat time.
        
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
        # d2: the d3 dataset with out the missing values 
        
        d1=dategen(k,M,varx=varx,icc=icc,iccm=iccm,intercept=intercept)
        d3=d1
        d3$y=ifelse(d3$r==1,NA,d3$y)
        d3$missing=d3$r
        d2=na.omit(d3)
        
        # calculate the weights for IPW
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
        
        miss=c(miss,sum(d1$r))  # save the missing percentages
        
        
        # Get the true effects.
        # the true effect of intervention arm is calculated by: 
        # apply GEE function on the full dataset without missing values (which is d1)
        # for 1000 times and the mean value of the 1000 results is calculated and is 
        # taken as the true margnial effects of the intervention arm.  
        
        ### True effect
        trues_ind=myTryCatch(geeglm(formula=y~arm,id=cluster, data = d1,
                                    family =  binomial("logit"),
                                    corstr = "independence"))
        trues_ex=myTryCatch(geeglm(formula=y~arm,id=cluster, data = d1,
                                   family =  binomial("logit"),
                                   corstr = "exchangeable"))
        
        ### Unadjusted CRA
        ucra_ind=myTryCatch(geeglm(formula=y~arm,id=cluster, data = d2,
                                   family =  binomial("logit"),
                                   corstr = "independence"))
        ucra_ex=myTryCatch(geeglm(formula=y~arm,id=cluster, data = d2,
                                  family =  binomial("logit"),
                                  corstr = "exchangeable"))
        
        ### Adjusted CRA
        cra_ind=myTryCatch(geeglm(formula=y~x+arm,id=cluster, data = d2,
                                  family =  binomial("logit"),
                                  corstr = "independence"))
        cra_ex=myTryCatch(geeglm(formula=y~x+arm,id=cluster, data = d2,
                                 family =  binomial("logit"),
                                 corstr = "exchangeable"))
        
        ### IPW without cluster effects 
        ipw_ind=myTryCatch(geeglm(formula=y~arm,id=cluster, data = d3,
                                  family =  binomial("logit"),
                                  weights = d3$weight,
                                  corstr = "independence"))
        ipw_ex=myTryCatch(geeglm(formula=y~arm,id=cluster, data = d3,
                                 family =  binomial("logit"),
                                 weights = d3$weight,
                                 corstr = "exchangeable"))
        
        ### IPW with cluster effects
        ipw_clu_ind=myTryCatch(geeglm(formula=y~arm,id=cluster, data = d3,
                                      family =  binomial("logit"),
                                      weights = d3$weight2,
                                      corstr = "independence"))
        ipw_clu_ex=myTryCatch(geeglm(formula=y~arm,id=cluster, data = d3,
                                     family =  binomial("logit"),
                                     weights = d3$weight2,
                                     corstr = "exchangeable"))
        
        
        
        ## Save results:
        # The myTryCatch function is used to avoid errors, which may stop the 
        # simulation loop
        # Also, the erros and warnings are stored in the results
        
        # the true effects with indenpent working correlation matrix
        if(is.null(trues_ind$value)==0){
          t1=summary(trues_ind$value)$coefficients['arm','Estimate'] # save the mean
          t2=summary(trues_ind$value)$coefficients['arm','Std.err']  # save the sd
          true_est_ind=c(true_est_ind,t1)
          true_std_ind=c(true_std_ind,t2)
          if(summary(trues_ind$value)$error==1){true_warn_ind=c(true_warn_ind,times)}
          if(summary(trues_ind$value)$error==0){true_warn_ind=c(true_warn_ind,0)}
        }
        if(is.null(trues_ind$value)==1){
          true_est_ind=c(true_est_ind,NA)
          true_std_ind=c(true_std_ind,NA)
          true_warn_ind=c(true_warn_ind,times)
        }
        
        # the true effects with exchangeable working correlation matrix
        if(is.null(trues_ex$value)==0){
          t1=summary(trues_ex$value)$coefficients['arm','Estimate']
          t2=summary(trues_ex$value)$coefficients['arm','Std.err']
          true_est_ex=c(true_est_ex,t1)
          true_std_ex=c(true_std_ex,t2)
          if(summary(trues_ex$value)$error==1){true_warn_ex=c(true_warn_ex,times)}
          if(summary(trues_ex$value)$error==0){true_warn_ex=c(true_warn_ex,0)}
        }
        if(is.null(trues_ex$value)==1){
          true_est_ex=c(true_est_ex,NA)
          true_std_ex=c(true_std_ex,NA)
          true_warn_ex=c(true_warn_ex,times)
        }
        
        # the unadjusted CRA with indenpent working correlation matrix
        if(is.null(ucra_ind$value)==0){
          t1=summary(ucra_ind$value)$coefficients['arm','Estimate']
          t2=summary(ucra_ind$value)$coefficients['arm','Std.err']
          ucra_est_ind=c(ucra_est_ind,t1)
          ucra_std_ind=c(ucra_std_ind,t2)
          if(summary(ucra_ind$value)$error==1){ucra_warn_ind=c(ucra_warn_ind,times)}
          if(summary(ucra_ind$value)$error==0){ucra_warn_ind=c(ucra_warn_ind,0)}
        }
        if(is.null(ucra_ind$value)==1){
          ucra_est_ind=c(ucra_est_ind,NA)
          ucra_std_ind=c(ucra_std_ind,NA)
          ucra_warn_ind=c(ucra_warn_ind,times)
        }
        
        # the adjusted CRA with indenpent working correlation matrix
        if(is.null(ucra_ex$value)==0){
          t1=summary(ucra_ex$value)$coefficients['arm','Estimate']
          t2=summary(ucra_ex$value)$coefficients['arm','Std.err']
          ucra_est_ex=c(ucra_est_ex,t1)
          ucra_std_ex=c(ucra_std_ex,t2)
          if(summary(ucra_ex$value)$error==1){ucra_warn_ex=c(ucra_warn_ex,times)}
          if(summary(ucra_ex$value)$error==0){ucra_warn_ex=c(ucra_warn_ex,0)}
        }
        if(is.null(ucra_ex$value)==1){
          ucra_est_ex=c(ucra_est_ex,NA)
          ucra_std_ex=c(ucra_std_ex,NA)
          ucra_warn_ex=c(ucra_warn_ex,times)
        }
        
        # the adjusted CRA with independent working correlation matrix
        if(is.null(cra_ind$value)==0){
          t1=summary(cra_ind$value)$coefficients['arm','Estimate']
          t2=summary(cra_ind$value)$coefficients['arm','Std.err']
          cra_est_ind=c(cra_est_ind,t1)
          cra_std_ind=c(cra_std_ind,t2)
          if(summary(cra_ind$value)$error==1){cra_warn_ind=c(cra_warn_ind,times)}
          if(summary(cra_ind$value)$error==0){cra_warn_ind=c(cra_warn_ind,0)}
        }
        if(is.null(cra_ind$value)==1){
          cra_est_ind=c(cra_est_ind,NA)
          cra_std_ind=c(cra_std_ind,NA)
          cra_warn_ind=c(cra_warn_ind,times)
        }
        
        #  the adjusted CRA with exchangeable working correlation matrix
        if(is.null(cra_ex$value)==0){
          t1=summary(cra_ex$value)$coefficients['arm','Estimate']
          t2=summary(cra_ex$value)$coefficients['arm','Std.err']
          cra_est_ex=c(cra_est_ex,t1)
          cra_std_ex=c(cra_std_ex,t2)
          if(summary(cra_ex$value)$error==1){cra_warn_ex=c(cra_warn_ex,times)}
          if(summary(cra_ex$value)$error==0){cra_warn_ex=c(cra_warn_ex,0)}
        }
        if(is.null(cra_ex$value)==1){
          cra_est_ex=c(cra_est_ex,NA)
          cra_std_ex=c(cra_std_ex,NA)
          cra_warn_ex=c(cra_warn_ex,times)
        }
        
        # the IPW with indendent working correlation matrix
        if(is.null(ipw_ind$value)==0){
          t1=summary(ipw_ind$value)$coefficients['arm','Estimate']
          t2=summary(ipw_ind$value)$coefficients['arm','Std.err']
          ipw_est_ind=c(ipw_est_ind,t1)
          ipw_std_ind=c(ipw_std_ind,t2)
          if(summary(ipw_ind$value)$error==1){ipw_warn_ind=c(ipw_warn_ind,times)}
          if(summary(ipw_ind$value)$error==0){ipw_warn_ind=c(ipw_warn_ind,0)}
        }
        if(is.null(ipw_ind$value)==1){
          ipw_est_ind=c(ipw_est_ind,NA)
          ipw_std_ind=c(ipw_std_ind,NA)
          ipw_warn_ind=c(ipw_warn_ind,times)
        }
        
        # the IPW with exchangeable working correlation matrix
        if(is.null(ipw_ex$value)==0){
          t1=summary(ipw_ex$value)$coefficients['arm','Estimate']
          t2=summary(ipw_ex$value)$coefficients['arm','Std.err']
          ipw_est_ex=c(ipw_est_ex,t1)
          ipw_std_ex=c(ipw_std_ex,t2)
          if(summary(ipw_ex$value)$error==1){ipw_warn_ex=c(ipw_warn_ex,times)}
          if(summary(ipw_ex$value)$error==0){ipw_warn_ex=c(ipw_warn_ex,0)}
        }
        if(is.null(ipw_ex$value)==1){
          ipw_est_ex=c(ipw_est_ex,NA)
          ipw_std_ex=c(ipw_std_ex,NA)
          ipw_warn_ex=c(ipw_warn_ex,times)
        }
        
        # the IPW with cluster effects with independent working correlation matrix
        if(is.null(ipw_clu_ind$value)==0){
          t1=summary(ipw_clu_ind$value)$coefficients['arm','Estimate']
          t2=summary(ipw_clu_ind$value)$coefficients['arm','Std.err']
          ipw_clu_est_ind=c(ipw_clu_est_ind,t1)
          ipw_clu_std_ind=c(ipw_clu_std_ind,t2)
          if(summary(ipw_clu_ind$value)$error==1){ipw_clu_warn_ind=c(ipw_clu_warn_ind,times)}
          if(summary(ipw_clu_ind$value)$error==0){ipw_clu_warn_ind=c(ipw_clu_warn_ind,0)}
        }
        if(is.null(ipw_clu_ind$value)==1){
          ipw_clu_est_ind=c(ipw_clu_est_ind,NA)
          ipw_clu_std_ind=c(ipw_clu_std_ind,NA)
          ipw_clu_warn_ind=c(ipw_clu_warn_ind,times)
        }
        
        # the IPW with cluster effects with exchangeable working correlation matrix
        if(is.null(ipw_clu_ex$value)==0){
          t1=summary(ipw_clu_ex$value)$coefficients['arm','Estimate']
          t2=summary(ipw_clu_ex$value)$coefficients['arm','Std.err']
          ipw_clu_est_ex=c(ipw_clu_est_ex,t1)
          ipw_clu_std_ex=c(ipw_clu_std_ex,t2)
          if(summary(ipw_clu_ex$value)$error==1){ipw_clu_warn_ex=c(ipw_clu_warn_ex,times)}
          if(summary(ipw_clu_ex$value)$error==0){ipw_clu_warn_ex=c(ipw_clu_warn_ex,0)}
        }
        if(is.null(ipw_clu_ex$value)==1){
          ipw_clu_est_ex=c(ipw_clu_est_ex,NA)
          ipw_clu_std_ex=c(ipw_clu_std_ex,NA)
          ipw_clu_warn_ex=c(ipw_clu_warn_ex,times)
        }
      }
      
      
      ## save the data in the data.frame
      
      # the estimated value of true effects
      est_ind=data.frame(true_est_ind=true_est_ind,ucra_est_ind=ucra_est_ind,cra_est_ind=cra_est_ind,
                         ipw_est_ind=ipw_est_ind,ipw_clu_est_ind=ipw_clu_est_ind)
      est_ex=data.frame(true_est_ex=true_est_ex,ucra_est_ex=ucra_est_ex,cra_est_ex=cra_est_ex,
                        ipw_est_ex=ipw_est_ex,ipw_clu_est_ex=ipw_clu_est_ex)
      
      # sd for the true effects
      std_ind=data.frame(true_std_ind=true_std_ind,ucra_std_ind=ucra_std_ind,cra_std_ind=cra_std_ind,
                         ipw_std_ind=ipw_std_ind,ipw_clu_std_ind=ipw_clu_std_ind)
      std_ex=data.frame(true_std_ex=true_std_ex,ucra_std_ex=ucra_std_ex,cra_std_ex=cra_std_ex,
                        ipw_std_ex=ipw_std_ex,ipw_clu_std_ex=ipw_clu_std_ex)
      
      # the warning times 
      warn_ind=data.frame(true_warn_ind=true_warn_ind,ucra_warn_ind=ucra_warn_ind,cra_warn_ind=cra_warn_ind,
                          ipw_warn_ind=ipw_warn_ind,ipw_clu_warn_ind=ipw_clu_warn_ind)
      warn_ex=data.frame(true_warn_ex=true_warn_ex,ucra_warn_ex=ucra_warn_ex,cra_warn_ex=cra_warn_ex,
                         ipw_warn_ex=ipw_warn_ex,ipw_clu_warn_ex=ipw_clu_warn_ex)
      
      # combine all the results and save in a RData file
      result=list(est_ind=est_ind,est_ex=est_ex,std_ind=std_ind,
                  std_ex=std_ex,warn_ind=warn_ind,warn_ex=warn_ex,miss=miss,maxw1=maxw1,maxw2=maxw2)
      names=paste('wrapk',k,icc,iccm,'.RData',sep='')
      save(result,file=names)
    }
  }
}
e1=Sys.time()
