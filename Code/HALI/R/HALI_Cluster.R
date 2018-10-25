#### HALI Data Analysis

### calculate the statistics for number of missingness in each cluster 

library(readstata13) # library for read in the dataset
library(jomo) # library for MMI
library(lme4) # library for GLM
library(geeM) # library for GEE

##### Read in dataset
setwd('/Users/yaolanqiu/Documents/HALI/DATA')
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

dat11 = dat1
dat11$missing = is.na(dat1$gll21_total)

dat_int = dat11[dat11$LIT_grp=='yes',]
dat_col = dat11[dat11$LIT_grp=='no',]

mean(aggregate(dat_int[, 'missing'], list(dat_int$school_id), sum)$x)
sd(aggregate(dat_int[, 'missing'], list(dat_int$school_id), sum)$x)

mean(aggregate(dat_col[, 'missing'], list(dat_col$school_id), sum)$x)
sd(aggregate(dat_col[, 'missing'], list(dat_col$school_id), sum)$x)

dat_int = dat0[dat0$LIT_grp=='yes',]
dat_col = dat0[dat0$LIT_grp=='no',]

pdf('Histogram of number of missing in each cluster.pdf')
par(mfrow=c(1,2))
hist(aggregate(dat_int[, 'missing'], list(dat_int$school_id), sum)$x,
     breaks =10, main = 'Intervention Arm',
     xlab = 'Number of missing in each cluster')
hist(aggregate(dat_col[, 'missing'], list(dat_col$school_id), sum)$x,breaks =10, 
     main = 'Control Arm',
     xlab = 'Number of missing in each cluster')
title(main="Histogram of number of missing in each cluster",outer=T)
dev.off()