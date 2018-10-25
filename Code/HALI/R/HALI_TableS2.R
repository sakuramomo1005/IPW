# HALI Data Analysis

#### Predictors of outcome (high literacy) and of missingness (of high literacy outcome) from HALI motivating data set.

# Draw Table S2

library(readstata13) # library for read in the dataset
library(lme4) # library for GLM
library(geeM)
library(jomo)

setwd("/Users/yaolanqiu/Documents/HALI/DATA")
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

## Predictors of outcome (high literacy) 
## Predictors of outcome missingness

#### model 1
# with outcome:
#  * high literacy

# with covariates:
#  * baseline score
#* age 
#* sex
#* with or without household education level
#* SES

model1_edu=glmer(y~LIT_grp+BL_gll21_total+age_child+sex
                 +schlevel_comp+BL_ses+(1|school_id),
                 data=dat1,
                 family=binomial('logit'),
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))

model1_edu2=glmer(y~LIT_grp+BL_gll21_total+age_child+sex
                  +BL_ses+(1|school_id),
                  data=dat1,
                  family=binomial('logit'),
                  control=glmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e5)))


#### Calculate the overall p value for household education level

edu=anova(model1_edu, model1_edu2, test="LRT")
round(edu$`Pr(>Chisq)`[2],3)

#### Similarly, calculate the overall p value for SES
model2_ses=glmer(y~LIT_grp+BL_gll21_total+age_child+sex+
                   schlevel_comp+BL_ses+(1|school_id),
                 data=dat1,
                 family=binomial('logit'),
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))
model2_ses2=glmer(y~LIT_grp+BL_gll21_total+age_child+sex+
                    schlevel_comp+(1|school_id),
                  data=dat1,
                  family=binomial('logit'),
                  control=glmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e5)))
ses=anova(model2_ses, model2_ses2, test="LRT")

#### Model 3 

# * outcome: missing indicator

model3_edu=glmer(missing~LIT_grp+BL_gll21_total+age_child+sex+schlevel_comp
                 +BL_ses+(1|school_id),
                 data=dat1,
                 family=binomial('logit'),
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))
model3_edu2=glmer(missing~LIT_grp+BL_gll21_total+age_child+sex
                  +BL_ses+(1|school_id),
                  data=dat1,
                  family=binomial('logit'),
                  control=glmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e5)))

# the overall p value for household education level
edu2=anova(model3_edu, model3_edu2, test="LRT")

#### Similarly, calculate the overall p value for SES
model4_ses=glmer(missing~LIT_grp+BL_gll21_total+age_child+sex+schlevel_comp
                 +BL_ses+(1|school_id),
                 data=dat1,
                 family=binomial('logit'),
                 control=glmerControl(optimizer="bobyqa",
                                      optCtrl=list(maxfun=2e5)))

model4_ses2=glmer(missing~LIT_grp+BL_gll21_total+age_child+sex+schlevel_comp
                  +(1|school_id),
                  data=dat1,
                  family=binomial('logit'),
                  control=glmerControl(optimizer="bobyqa",
                                       optCtrl=list(maxfun=2e5)))
ses2=anova(model4_ses, model4_ses2, test="LRT")


#### Draw Table S2
table1=round(summary(model2_ses)$coefficient,3)
table2=round(summary(model4_ses)$coefficient,3)
table1=cbind(formatC(table1[,c(1:3)],2,format="f"),table1[,4])
table1[,4][table1[,4]=="0"] = '<0.001'

table2=cbind(formatC(table2[,c(1:3)],2,format="f"),table2[,4])
table2[,4][table2[,4]=="0"] = '<0.001'

ass1=rbind(table1[1:5,],
           c('','','',formatC(round(edu$`Pr(>Chisq)`[2],3),3,format="f")),
           table1[6:8,],
           c("","","",formatC(round(ses$`Pr(>Chisq)`[2],3),3,format="f")),
           table1[9:12,]
)

ass2=rbind(table2[1:5,],
           c('','','',formatC(round(edu2$`Pr(>Chisq)`[2],3),3,format="f")),
           table2[6:8,],
           c("","","",formatC(round(ses2$`Pr(>Chisq)`[2],3),3,format="f")),
           table2[9:12,]
) 
ass=cbind(c(
  'Intercept','Intervention','Baseline literacy score',
  'Age','Sex (Female)','HH education','  primary','  secondary',
  '  college/degree','SES household','  Poor','  Median Poor','  Less Poor',
  '  Least Poor'
),ass1,ass2)
rownames(ass)=NULL
ass = ass[,c(1,2,3,5,6,7,9)]
ass[c(7:9,11:14),c(4,7)] = ''
colnames(ass) = c('','Estimate','SE','p-value','Estimate','SE','p-value')

tableS2 = ass
tableS2