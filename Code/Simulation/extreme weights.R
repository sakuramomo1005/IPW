
#install.packages('lme4')
#install.packages('geepack')

library(lme4)
library(geepack)


#####################################################################################
#                                    Functions
#####################################################################################


# function to catch errors and warns
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

# function to calculate the variance based on ICC
missing_icc=function(icc){
  pi=3.142
  a=pi^2/(3*(1/icc-1))
  return(a)
}

# function to calculate the missingness percentage in dataset
missing_per=function(data){
  res=sum(data$r)/dim(data)[1]
  return(res)
}

# expit function
expit=function(x){y=exp(x)/(1+exp(x));return(y)}

# data generation function
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
##                                 Count the weights
#####################################################################################

##### Calculate how many weights exceed 20, 100, 500, and 1000

# I saved the weigths and then calculated the 
# percentages in case we may need more weigths characteristics 
## Extreme weigths calculation

# calculate the weigths within simulations for W-GEE and CW-GEE

T=1000;varx=0.2;M=50

s1=Sys.time()

for(k in c(25,50)){
  for(icc in c(0.05,0.1,0.2)){
    for(iccm in c(0,0.1,0.3,0.5)){
      
      weight=c()
      weight2=c()
      ts=c() # save the times 
      
      for(times in 1:T){
        print(paste('k',k,'icc',icc,'iccm',iccm,'times',times))
        set.seed(times)
        
        #tune the intercept
        mis=c()
        for(intercept in seq(-5,5,0.1)){
          temp=dategen(k,M,varx=varx,icc=icc,iccm=iccm,intercept=intercept)  
          mis=c(mis,missing_per(temp))
        }
        
        intercept=seq(-5,5,0.1)[which.min(abs(mis-0.3))]
        
        #
        d1=dategen(k,M,varx=varx,icc=icc,iccm=iccm,intercept=intercept)
        d3=d1
        d3$y=ifelse(d3$r==1,NA,d3$y)
        d3$missing=d3$r
        
        w1=glm(missing ~x + arm, data = d3,
               family = binomial(link='logit'))
        w2=glmer(missing ~x + arm+(1|cluster) , data = d3,
                 family = binomial(link='logit'))
        
        w1=predict(w1,type="response")
        w2=predict(w2,type="response")
        
        weight=c(weight,1/w1)
        weight2=c(weight2,1/w2)
        ts=c(ts,rep(times,length(1/w1)))
      }
      result=data.frame(ts=ts,weight=weight,weight2=weight2)
      names=paste('wkicciccm',k,icc,iccm,'.RData',sep='')
      save(result,file = names)
    }
  }
}

s2=Sys.time()

