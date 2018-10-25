##### HALI Data Analysis
##### Five methods

##### Draw Table 1b

#Estimated intervention effects 
#under five different GEE approaches 
#(with robust SE and exchangeable working correlation matrix)

library(readstata13) # library for read in the dataset
library(lme4) # library for GLM
library(geeM)
library(jomo)

### The pool function for MMI
mypool <- function(mean0,sd0,num=5,J=50){
  # input: 
  # mean0: a vector of the values of estimated beta parameter, with length equals to the imputation time.
  # sd0: a vector of the standard deviations of the estimated beta parameter, with length equals to the imputation time.
  # num: impuation time
  # J: the number of clusters in each intervention arm 
  
  # count the times of NA in the input vector.
  na_times <- sum(is.na(mean0)) 
  # the number of imputations without NA. 
  num_actual <- num-na_times 
  
  # the MI estimate of the beta parameter
  m <- mean(mean0,na.rm=TRUE) 
  
  # estimate of average wihtin-imputation variance 
  # i.e. based on the SE^2 of the beta parameter from each fitted model
  W <- mean(sd0^2,na.rm=TRUE) 
  
  # estimate of between-imputation variance 
  # i.e. empirical SD of the point estimates of beta parameter
  B <- var(mean0,na.rm=TRUE) 
  
  # estimate of total variance 
  # i.e. will need to take the sqrt of this to use for inference in making confidence intervals etc.
  v_hat <- W+(1+1/num_actual)*B 
  v_hat_sqrt<-sqrt(v_hat)
  
  # Confindence interval
  # Testing based on naive normal distribution assumption
  l <- m-1.96*v_hat_sqrt
  u <- m+1.96*v_hat_sqrt
  
  # Testing based on standard results from MI literature
  # i.e. df of t distribution for testing based on standard results from MI literature
  df_denom <- (1+1/num_actual)*B
  df_part <- 1+W/df_denom
  df_t <- (num_actual-1)*df_part^2 # adjusted df of t distribution
  l_t <- m-qt(0.975,df_t)*v_hat_sqrt
  u_t <- m+qt(0.975,df_t)*v_hat_sqrt
  
  # Testing based on results from MMI literature, Barnard and Rubin (1999), 
  # df of t distribution for testing based on results in re adjustment for MMI feature
  
  #df for full data
  df_com <- 2*J - 2 
  
  # calculate the adjusted df based on the literature.
  parenthesis <- 1+df_denom*(1/W) 
  df_obs <- df_com*((df_com+1)/(df_com+3))*(1/parenthesis) 
  df_adj_t <- 1/(1/df_t + 1/df_obs) 
  
  # Print the results
  print("Standard t df and Barnard/Rubin adjusted t df");
  print(c(df_t, df_adj_t))
  print("97.5% quantiles from standard t df and Barnard/Rubin adjusted t df");
  print(c(qt(0.975,df_t), qt(0.975,df_adj_t)))
  
  # The confidence interval calculated based on df of standard t distribution and adjusted df of t distribution 
  l_t <- m-qt(0.975,df_t)*v_hat_sqrt
  u_t <- m+qt(0.975,df_t)*v_hat_sqrt
  l_adj_t <- m-qt(0.975,df_adj_t)*v_hat_sqrt
  u_adj_t <- m+qt(0.975,df_adj_t)*v_hat_sqrt
  
  return(list(mean=m,std=v_hat_sqrt,
              df_t=df_t,
              df_adj_t=df_adj_t))
}

### Read in and clean data
setwd(' ')
dat = read.dta13('HALI_CLASS1_2539_MASTER_EDUC_FU1_FU2_LONG_accounts_withdrawals_BL_AS_COV_AND_OTHER_COV_17.11.2016.dta',
                 nonint.factors = TRUE, generate.factors=TRUE)
# the full data named dat

### We only focused on 9 month data
data0 <- dat[dat$visit=='9-month FU',]
dim(data0) # 2539 220
# the data with only 9 month named data0

### Select the input variables
dat0 <- data0[,c("school_id","LIT_grp",
                 "BL_gll21_total", 'age_child','sex',
                 'schlevel_comp',"BL_ses",
                 "gll21_total")]
rownames(dat0) <- NULL
dim(dat0) # 2539 8

### Remove the data with missing covariates 
dat1 = dat0[(is.na(dat0$school_id) == 0 & 
               is.na(dat0$LIT_grp) == 0 &
               is.na(dat0$BL_gll21_total) == 0 &
               is.na(dat0$age_child)== 0 & 
               is.na(dat0$sex) == 0 & 
               is.na(dat0$schlevel_comp) == 0 &
               is.na(dat0$BL_ses) ==0) ,]
dim(dat1) #  2465 8

# Add misisng indicator
dat1$missing = is.na(dat1$gll21_total)
# Add high literacy indicator
y = ifelse(dat1$gll21_total>10,1,0)
dat1$y = y

#### draw table 1 b

# formulas
# formula for GEE
formula1 = y ~ LIT_grp
formula2 = y ~ LIT_grp + BL_gll21_total + 
  age_child + sex + schlevel_comp + BL_ses

# formula to calculate the weights 
formula3 = missing ~ LIT_grp + BL_gll21_total + 
  age_child + sex + schlevel_comp + 
  BL_ses
formula4 = missing ~ LIT_grp + BL_gll21_total + 
  age_child + sex + schlevel_comp + 
  BL_ses + (1 | school_id)


# Empty vectors that save results later.
means <- c();sds <- c();cis <- c();ps <- c()

################ CRA-GEE ################ 
# independent working correlation matrix
ucra_ind <- geem(formula1,id = school_id, 
                 data = dat1,
                 family =  binomial("logit"),
                 corstr = "independence")
summary(ucra_ind)

# save the results 
m <- summary(ucra_ind)$beta[2]
s <- summary(ucra_ind)$se.robust[2]
p <- summary(ucra_ind)$p[2]
ciu <- m-1.96*s; cil <- m+1.96*s

# save results 
means <- c(means,exp(m))
sds <- c(sds,s)
ps <- c(ps,p)
ciu <- formatC(exp(ciu),2,format="f"); cil <- formatC(exp(cil),2,format="f")
cis <- c(cis,paste('(',ciu,', ',cil,')',sep=''))

# exchangeable working correlation matrix
ucra_ex <- geem(formula = formula1,id = school_id, 
                data = dat1,
                family =  binomial("logit"),
                corstr = "exchangeable")
ucra_ex
summary(ucra_ex)
# save the results 
m <- summary(ucra_ex)$beta[2]
s <- summary(ucra_ex)$se.robust[2]
p <- summary(ucra_ex)$p[2]
ciu <- m-1.96*s; cil <- m+1.96*s

# save results 
means <- c(means,exp(m))
sds <- c(sds,s)
ps <- c(ps,p)
ciu <- formatC(exp(ciu),2,format="f"); cil <- formatC(exp(cil),2,format="f")
cis <- c(cis,paste('(',ciu,', ',cil,')',sep=''))

################ A-CRA-GEE ################ 
# independent working correlation matrix
cra_ind <- geem(formula = formula2,
                id = school_id, data = dat1,
                family =  binomial("logit"),
                corstr = "independence")
summary(cra_ind)
# save the results 
m <- summary(cra_ind)$beta[2]
s <- summary(cra_ind)$se.robust[2]
p <- summary(cra_ind)$p[2]
ciu <- m-1.96*s; cil <- m+1.96*s

# save results 
means <- c(means,exp(m))
sds <- c(sds,s)
ps <- c(ps,p)
ciu <- formatC(exp(ciu),2,format="f"); cil <- formatC(exp(cil),2,format="f")
cis <- c(cis,paste('(',ciu,', ',cil,')',sep=''))

# exchangeable working correlation matrix
cra_ex <- geem(formula = formula2,
               id = school_id, data = dat1,
               family =  binomial("logit"),
               corstr = "exchangeable")
cra_ex
summary(cra_ex)
# save the results 
m <- summary(cra_ex)$beta[2]
s <- summary(cra_ex)$se.robust[2]
p <- summary(cra_ex)$p[2]
ciu <- m-1.96*s; cil <- m+1.96*s

# save results 
means <- c(means,exp(m))
sds <- c(sds,s)
ps <- c(ps,p)
ciu <- formatC(exp(ciu),2,format="f"); cil <- formatC(exp(cil),2,format="f")
cis <- c(cis,paste('(',ciu,', ',cil,')',sep=''))

################ W-GEE ################ 

# 1. calculate the weights
w1 <- glm(formula = formula3, data = dat1,
          family = binomial(link='logit'))
w2 <- glmer(formula = formula4, data = dat1,
            family = binomial(link='logit'),
            control=glmerControl(optimizer="bobyqa",
                                 optCtrl=list(maxfun=2e5)))
# increase the iteration times to avoid non-convergence

w1 <- predict(w1,type="response")  # get the weights value from the glm
w2 <- predict(w2,type="response")  # get the weights value form the glmer

# add the weights to the dataset
dat1$weight <- 1/w1
dat1$weight2 <- 1/w2

# with independent working correlation matrix
ipw_ind <- geem(formula = y~LIT_grp,
                id = school_id, data = dat1,
                family =  binomial("logit"),
                weights = dat1$weight,
                corstr = "independence")
summary(ipw_ind)

# save the results 
m <- summary(ipw_ind)$beta[2]
s <- summary(ipw_ind)$se.robust[2]
p <- summary(ipw_ind)$p[2]
ciu <- m-1.96*s; cil <- m+1.96*s

# save results 
means <- c(means,exp(m))
sds <- c(sds,s)
ps <- c(ps,p)
ciu <- formatC(exp(ciu),2,format="f"); cil <- formatC(exp(cil),2,format="f")
cis <- c(cis,paste('(',ciu,', ',cil,')',sep=''))

# with exchangeable working correlation matrix
ipw_ex <- geem(formula = y~LIT_grp,
               id = school_id, data = dat1,
               family =  binomial("logit"),
               weights = dat1$weight,
               corstr = "exchangeable")
summary(ipw_ex)

# save the results 
m <- summary(ipw_ex)$beta[2]
s <- summary(ipw_ex)$se.robust[2]
p <- summary(ipw_ex)$p[2]
ciu <- m-1.96*s; cil <- m+1.96*s

# save results 
means <- c(means,exp(m))
sds <- c(sds,s)
ps <- c(ps,p)
ciu <- formatC(exp(ciu),2,format="f"); cil <- formatC(exp(cil),2,format="f")
cis <- c(cis,paste('(',ciu,', ',cil,')',sep=''))


################ CW-GEE ################
# independent working correlation matrix
ipw_clu_ind <- geem(formula = y~LIT_grp,
                    id=school_id, data = dat1,
                    family =  binomial("logit"),
                    weights = dat1$weight2,
                    corstr = "independence")
summary(ipw_clu_ind)

# save the results 
m <- summary(ipw_clu_ind)$beta[2]
s <- summary(ipw_clu_ind)$se.robust[2]
p <- summary(ipw_clu_ind)$p[2]
ciu <- m-1.96*s; cil <- m+1.96*s
# save the results 
means <- c(means,exp(m))
ps <- c(ps,p)
sds <- c(sds,s)
ciu = exp(ciu); cil = exp(cil)
cis <- c(cis,paste('(',round(ciu,2),', ',round(cil,2),')',sep=''))

# exchangeable working correlation matrix
ipw_clu_ex <- geem(formula = y~LIT_grp,
                   id = school_id, data = dat1,
                   family =  binomial("logit"),
                   weights = dat1$weight2,
                   corstr = "exchangeable")
summary(ipw_clu_ex)

# save the results 
m <- summary(ipw_clu_ex)$beta[2]
s <- summary(ipw_clu_ex)$se.robust[2]
p <- summary(ipw_clu_ex)$p[2]
ciu <- m-1.96*s; cil <- m+1.96*s

# save results 
means <- c(means,exp(m))
sds <- c(sds,s)
ps <- c(ps,p)
ciu <- formatC(exp(ciu),2,format="f"); cil <- formatC(exp(cil),2,format="f")
cis <- c(cis,paste('(',ciu,', ',cil,')',sep=''))

################ MMI-GEE ################ 
covariates <- c("BL_gll21_total","age_child","sex","schlevel_comp","BL_ses")
data.miss <- dat1
# data frame for response variables with missing values
y.cat <- data.frame(outcome=factor(data.miss$y))  
# number of levels in outcome variable
y.numcat <- c(2)     
# data frame for clusters
clus <- data.frame(clus=data.miss$school_id)         
# covariates, includes 1s.
nobs <- dim(data.miss)[1]
x <- data.frame(rep(1,nobs),
                data.miss[,covariates],
                group=data.miss$LIT_grp)

# run to generate Nimp full datasets
# imp is full datasets 
# Nimp <- 15
# generate the complete dataset

set.seed(123)
Nimp <- 15
imp <- jomo1rancat(Y.cat = y.cat, Y.numcat = y.numcat, X = x,
                   clus = clus, nburn = 100, nbetween = 25, 
                   nimp = Nimp,output = 0)

mmi_est_ind <- c();mmi_est_ex <- c()
mmi_std_ind <- c();mmi_std_ex <- c()
mmi_p_ind <- c();mmi_p_ex <- c()

for(i in 1:Nimp){
  temp <- imp[imp$Imputation==i,]
  rownames(temp) <- NULL
  temp$outcome <- as.numeric(temp$outcome)-1
  
  # run the analysis for Nimp times
  mmi_ind <- geem(formula=outcome~group,
                  id=clus , data = temp,
                  family =  binomial("logit"),
                  corstr = "independence")
  mmi_ex <- geem(formula=outcome~group,
                 id=clus , data = temp,
                 family =  binomial("logit"),
                 corstr = "exchangeable")
  
  # save the results 
  mmi_est_ind <- c(mmi_est_ind,
                   summary(mmi_ind)$beta[2])
  mmi_std_ind <- c(mmi_std_ind,
                   summary(mmi_ind)$se.robust[2])
  
  # save the results 
  mmi_est_ex <- c(mmi_est_ex,
                  summary(mmi_ex)$beta[2])
  mmi_std_ex <- c(mmi_std_ex,
                  summary(mmi_ex)$se.robust[2])
  mmi_p_ind <- c(mmi_p_ind,summary(mmi_ind)$p[2])
  mmi_p_ex <- c(mmi_p_ex,summary(mmi_ex)$p[2])
}

# pool results

pool1 <- mypool(mmi_est_ind,mmi_std_ind,num=Nimp)
pool2 <- mypool(mmi_est_ex,mmi_std_ex,num=Nimp)

# save the results 
# independent working correlation matrix
ciu <- pool1$mean - qt(0.975, pool1$df_adj_t)*pool1$std;
cil <- pool1$mean + qt(0.975, pool1$df_adj_t)*pool1$std
means <- c(means,exp(pool1$mean))
ps <- c(ps,round(mean(mmi_p_ind),3))
sds <- c(sds,pool1$std)
ciu <- formatC(exp(ciu),2,format="f"); cil <- formatC(exp(cil),2,format="f")
cis <- c(cis,paste('(',ciu,', ',cil,')',sep=''))

# exchangeable working correlation matrix
ciu <- pool2$mean - qt(0.975, pool2$df_adj_t)*pool2$std;
cil <- pool2$mean + qt(0.975, pool2$df_adj_t)*pool2$std
# save the results
means <- c(means,exp(pool2$mean))
sds <- c(sds,pool2$std) 
ps <- c(ps,round(mean(mmi_p_ex),3))
ciu <- formatC(exp(ciu),2,format="f"); cil <- formatC(exp(cil),2,format="f")
cis <- c(cis,paste('(',ciu,', ',cil,')',sep=''))



################ Draw table ################
# use the saved resutls to draw summary table
# use the saved resutls to draw summary table

coefs_ind <- formatC(means[seq(1,10,2)],2,format="f")
coefs_ex <- formatC(means[seq(2,10,2)],2,format="f")
std_ind <- formatC(sds[seq(1,10,2)],2,format="f")
std_ex <- formatC(sds[seq(2,10,2)],2,format="f")
cis_ind <- formatC(cis[seq(1,10,2)],2,format="f")
cis_ex <- formatC(cis[seq(2,10,2)],2,format="f")
p_ind <- formatC(ps[seq(1,10,2)],3,format="f")
p_ex <- formatC(ps[seq(2,10,2)],3,format="f")
names <- c('CRA-GEE','A-CRA-GEE','W-GEE','CW-GEE','MMI-GEE')

table_sum <- cbind(names,coefs_ex,std_ex,cis_ex,p_ex,coefs_ind,std_ind,cis_ind,p_ind)
colnames(table_sum) <- c('Method','Estimated','SE','Confidence Interval','p-value',
                         'Estimated','SE','Confidence Interval','p-value')
table_sum[table_sum[,c(5)] == "0.000",5] <- '<0.001'
table_sum[table_sum[,c(9)] == "0.000",9] <- '<0.001'
table_sum

