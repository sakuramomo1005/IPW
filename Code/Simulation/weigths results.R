# Draw results table for the counts of extreme weights 

w=c()

k = 50; icc = 0.2; iccm = 0.5

#### weights
for(k in c(10, 25, 50)){
  for(icc in c(0.05, 0.1, 0.2)){
    for(iccm in c(0, 0.1, 0.3, 0.5)){
      
      names=paste('wkicciccm',k,icc,iccm,'.RData',sep='')
      load(names)
      print(names)

      w1_20=c();w2_20=c()
      for(i in 1:1000){
        if(i %% 200==0){print(i)}
        temp=result[result$ts==i,]
        w1_20=c(w1_20,sum(temp$weight>20))
        w2_20=c(w2_20,sum(temp$weight2>20))
      }
      
      w1_100=c();w2_100=c()
      for(i in 1:1000){
        if(i %% 200==0){print(i)}
        temp=result[result$ts==i,]
        w1_100=c(w1_100,sum(temp$weight>100))
        w2_100=c(w2_100,sum(temp$weight2>100))
      }
      
      w1_500=c();w2_500=c()
      for(i in 1:1000){
        if(i %% 200==0){print(i)}
        temp=result[result$ts==i,]
        w1_500=c(w1_500,sum(temp$weight>500))
        w2_500=c(w2_500,sum(temp$weight2>500))
      }
      
      w1_1000=c();w2_1000=c()
      for(i in 1:1000){
        if(i %% 200==0){print(i)}
        temp=result[result$ts==i,]
        w1_1000=c(w1_1000,sum(temp$weight>1000))
        w2_1000=c(w2_1000,sum(temp$weight2>1000))
      }
      
      temp=c(mean(w1_20),mean(w1_100),mean(w1_500),mean(w1_1000),
             mean(w2_20),mean(w2_100),mean(w2_500),mean(w2_1000))
      w=rbind(w,c(k,icc,iccm,temp))
      print(w)
    }
  }
}

colnames(w) = c('k','icc','iccm','indw20','indw100','indw500','indw1000',
                'exw20','exw100','exw500','exw1000')
iccm = rep(c(0, 0.1, 0.3, 0.5),9)
icc = rep(rep(c(0.05, 0.1, 0.2),each = 4),3)
k = rep(c(10,25,50),each = 12)

W = cbind(k,icc,iccm,w)

W