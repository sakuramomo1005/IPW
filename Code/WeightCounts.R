#####################################################################################
#                                    Roadmap
# This is the R code for:
#  Get the summary of weights. Count how many of weights that exceed 20, 100, 500, 1000.

# Details:
# The simulation steps are:
#  1. Generate a full datasets without missingness
#  2. Generate missing values in outcome based on the function we mention in the paper
#  3. Calculate the weights of missingness
#  4. Save the results 
#  5. Repeat step 1-4 for 1000 times
#  6. Load in the saved RData files and counts how many weigths that exceed 20, 100,
#     500, and 1000. Get the mean value for each of these features in the 
#     1000 replicates
#####################################################################################

library(lme4)
library(geepack)
library(knitr)


T=1000;varx=0.2;M=50

s1=Sys.time()

for(k in c(10,25,50)){
  for(icc in c(0.05,0.1,0.2)){
    for(iccm in c(0,0.1,0.3,0.5,0.7,0.9)){
      
      weight=c() # the weights without cluster effects
      weight2=c() # the weight with cluster effects
      ts=c() # save the times of simulation
      
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



# Analyze the weights and make a table

w=c()
for(k in c(10, 25, 50)){
  for(icc in c(0.05, 0.1, 0.2)){
    for(iccm in c(0, 0.1, 0.3, 0.5, 0.7, 0.9)){
      
      # load in the saved data
      names=paste('wkicciccm',k,icc,iccm,'.RData',sep='')
      load(names)
      print(names)
      
      # weights larger than 20 in each replicate
      w1_20=c();w2_20=c()
      for(i in 1:1000){
        if(i %% 200==0){print(i)}
        temp=result[result$ts==i,]
        w1_20=c(w1_20,sum(temp$weight>20))
        w2_20=c(w2_20,sum(temp$weight2>20))
      }
      
      # weights larger than 100 in each replicate
      w1_100=c();w2_100=c()
      for(i in 1:1000){
        if(i %% 200==0){print(i)}
        temp=result[result$ts==i,]
        w1_100=c(w1_100,sum(temp$weight>100))
        w2_100=c(w2_100,sum(temp$weight2>100))
      }
      
      # weights larger than 500 in each replicate
      w1_500=c();w2_500=c()
      for(i in 1:1000){
        if(i %% 200==0){print(i)}
        temp=result[result$ts==i,]
        w1_500=c(w1_500,sum(temp$weight>500))
        w2_500=c(w2_500,sum(temp$weight2>500))
      }
      
      # weights larger than 1000 in each replicate
      w1_1000=c();w2_1000=c()
      for(i in 1:1000){
        if(i %% 200==0){print(i)}
        temp=result[result$ts==i,]
        w1_1000=c(w1_1000,sum(temp$weight>1000))
        w2_1000=c(w2_1000,sum(temp$weight2>1000))
      }
      
      # save the results
      temp=c(mean(w1_20),mean(w2_20),mean(w1_100),mean(w2_100),
             mean(w1_500),mean(w2_500),mean(w1_1000),mean(w2_1000))
      w=rbind(w,c(k,icc,iccm,temp))
      print(w)
    }
  }
}

# Draw the table: the percentage of weights 
ww=w[,4:11]
www=cbind(w[,1:3],round(ww/5000,2))
colnames(www)=c('K','ICC','ICCM','W1_20','W2_20','W1_100','W2_100','W1_500',
                'W2_500','W1_1000','W2_1000')
kable(www)

# Draw the table: the counts of weights
ww=w[,4:11]
www=cbind(w[,1:3],round(ww))
colnames(www)=c('K','ICC','ICCM','W1_20','W2_20','W1_100','W2_100','W1_500',
                'W2_500','W1_1000','W2_1000')
kable(www)
