# HALI Data Analysis

# Baseline and outcome characteristics 
# of motivating HALI data set for
# n=2465 participants with complete baseline covariates

#### Draw Table 1A

##### Library
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


##### Draw the table 1:

# Number of student in intervention arm and control arm:
dim(dat1[dat1$LIT_grp=='yes',]) # 1230
dim(dat1[dat1$LIT_grp=='no',])  # 1235

col1 = col2 =col3 =c()
### Baseline cluster characteristics:
col1 = c(col1, 'Baseline cluster characteristics')
col2 = c(col2, '')
col3 = c(col3, '')
## cluster number:
length(table(dat1[dat1$LIT_grp=='yes',]$school_id)) # intervention 51
length(table(dat1[dat1$LIT_grp=='no',]$school_id)) # control 50

col1 = c(col1,'Number')
col2 = c(col2, length(table(dat1[dat1$LIT_grp=='yes',]$school_id)) )
col3 = c(col3, length(table(dat1[dat1$LIT_grp=='no',]$school_id)))

## cluster size (mean (SD)):
# intervention 
# cell 24.1 (3.3)
paste(round(mean(table(dat1[dat1$LIT_grp=='yes',]$school_id)),1),' (',
      round(sd(table(dat1[dat1$LIT_grp=='yes',]$school_id)),1),")",sep='')
# control 
# cell 24.7 (1.9)
paste(round(mean(table(dat1[dat1$LIT_grp=='no',]$school_id)),1),' (',
      round(sd(table(dat1[dat1$LIT_grp=='no',]$school_id)),1),")",sep='')

col1 = c(col1,'Cluster Size – mean (SD)')
col2 = c(col2, paste(round(mean(table(dat1[dat1$LIT_grp=='yes',]$school_id)),1),' (',
                     round(sd(table(dat1[dat1$LIT_grp=='yes',]$school_id)),1),")",sep=''))
col3 = c(col3, paste(round(mean(table(dat1[dat1$LIT_grp=='no',]$school_id)),1),' (',
                     round(sd(table(dat1[dat1$LIT_grp=='no',]$school_id)),1),")",sep=''))


### Baseline child-level characteristics:
col1 = c(col1, 'Baseline child-level characteristics - % (n)')
col2 = c(col2, '')
col3 = c(col3, '')
# intervention
# gender - female: 47.9% (589)
# 589
dim(dat1[dat1$sex == 'Female' & dat1$LIT_grp == 'yes',])[1] 
# 47.9%
100 * round(dim(dat1[dat1$sex == 'Female' & 
                 dat1$LIT_grp == 'yes',])[1]/ dim(dat1[dat1$LIT_grp == 'yes',])[1],3)
# cell 47.9% (589)
paste(100 * round(dim(dat1[dat1$sex == 'Female' & 
                             dat1$LIT_grp == 'yes',])[1]/ dim(dat1[dat1$LIT_grp == 'yes',])[1],3),
      "% (",dim(dat1[dat1$sex == 'Female' & dat1$LIT_grp == 'yes',])[1],")",sep='')

# control
# gender - female:  49.5% (611)
# 611
dim(dat1[dat1$sex == 'Female' & dat1$LIT_grp == 'no',])[1] 
# 49.5%
100 * round(dim(dat1[dat1$sex == 'Female' & 
                       dat1$LIT_grp == 'no',])[1]/ dim(dat1[dat1$LIT_grp == 'no',])[1],3)
# cell 47.9% (589)
paste(100 * round(dim(dat1[dat1$sex == 'Female' & 
                             dat1$LIT_grp == 'no',])[1]/ dim(dat1[dat1$LIT_grp == 'no',])[1],3),
      "% (",dim(dat1[dat1$sex == 'Female' & dat1$LIT_grp == 'no',])[1],")",sep='')

col1 = c(col1, 'Female')
col2 = c(col2, paste(100 * round(dim(dat1[dat1$sex == 'Female' & 
                                            dat1$LIT_grp == 'yes',])[1]/ dim(dat1[dat1$LIT_grp == 'yes',])[1],3),
                     "% (",dim(dat1[dat1$sex == 'Female' & dat1$LIT_grp == 'yes',])[1],")",sep=''))
col3 = c(col3, paste(100 * round(dim(dat1[dat1$sex == 'Female' & 
                                            dat1$LIT_grp == 'no',])[1]/ dim(dat1[dat1$LIT_grp == 'no',])[1],3),
                     "% (",dim(dat1[dat1$sex == 'Female' & dat1$LIT_grp == 'no',])[1],")",sep=''))

# intervention
# age - mean(sd)
round(mean(dat1[dat1$LIT_grp == 'yes',]$age_child),1) # 7.1
round(sd(dat1[dat1$LIT_grp == 'yes',]$age_child),1) # 1.7
# cell 7.7 (1.7)
paste(round(mean(dat1[dat1$LIT_grp == 'yes',]$age_child),1),
      ' (', round(sd(dat1[dat1$LIT_grp == 'yes',]$age_child),1),
      ')', sep = '')

# control
# age - mean(sd)
round(mean(dat1[dat1$LIT_grp == 'no',]$age_child),1) # 7.9
round(sd(dat1[dat1$LIT_grp == 'no',]$age_child),1) # 1.7
# cell 7.9 (1.7)
paste(round(mean(dat1[dat1$LIT_grp == 'no',]$age_child),1),
      ' (', round(sd(dat1[dat1$LIT_grp == 'no',]$age_child),1),
      ')', sep = '')

col1 = c(col1, 'Age – mean (sd)')
col2 = c(col2, paste(round(mean(dat1[dat1$LIT_grp == 'yes',]$age_child),1),
                     ' (', round(sd(dat1[dat1$LIT_grp == 'yes',]$age_child),1),
                     ')', sep = ''))
col3 = c(col3, paste(round(mean(dat1[dat1$LIT_grp == 'no',]$age_child),1),
                     ' (', round(sd(dat1[dat1$LIT_grp == 'no',]$age_child),1),
                     ')', sep = ''))

# intervention
# household education:

col1 = c(col1, 'Household head education')
col2 = c(col2, '')
col3 = c(col3, '')

table(dat1[dat1$LIT_grp=='yes',]$schlevel_comp)
100 * round(table(dat1[dat1$LIT_grp=='yes',]$schlevel_comp)/dim(dat1[dat1$LIT_grp=='yes',])[1],3)
# "29.1% (358)" "55.6% (684)" "11.5% (141)" "3.8% (47)" 
paste(100 * round(table(dat1[dat1$LIT_grp=='yes',]$schlevel_comp)/dim(dat1[dat1$LIT_grp=='yes',])[1],3),
      "% (",
      table(dat1[dat1$LIT_grp=='yes',]$schlevel_comp), ")",sep='')

# control
# household education:
table(dat1[dat1$LIT_grp=='no',]$schlevel_comp)
100 * round(table(dat1[dat1$LIT_grp=='no',]$schlevel_comp)/dim(dat1[dat1$LIT_grp=='no',])[1],3)
# "34.4% (425)" "52.7% (651)" "10.6% (131)" "2.3% (28)" 
paste(100 * round(table(dat1[dat1$LIT_grp=='no',]$schlevel_comp)/dim(dat1[dat1$LIT_grp=='no',])[1],3),
      "% (",
      table(dat1[dat1$LIT_grp=='no',]$schlevel_comp), ")",sep='')

col1 = c(col1,c('Did not complete primary education',
                'Primary',
                'Secondary',
                'College/degree'))
col2 = c(col2, paste(100 * round(table(dat1[dat1$LIT_grp=='yes',]$schlevel_comp)/dim(dat1[dat1$LIT_grp=='yes',])[1],3),
                     "% (",
                     table(dat1[dat1$LIT_grp=='yes',]$schlevel_comp), ")",sep=''))
col3 = c(col3, paste(100 * round(table(dat1[dat1$LIT_grp=='no',]$schlevel_comp)/dim(dat1[dat1$LIT_grp=='no',])[1],3),
                     "% (",
                     table(dat1[dat1$LIT_grp=='no',]$schlevel_comp), ")",sep=''))

# intervention
# SES
col1 = c(col1, 'Household socioeconomic status (SES)')
col2 = c(col2, '')
col3 = c(col3, '')

table(dat1[dat1$LIT_grp=='yes',]$BL_ses)
100 * round(table(dat1[dat1$LIT_grp=='yes',]$BL_ses)/dim(dat1[dat1$LIT_grp=='yes',])[1],3)
# "19% (234)"   "19.9% (245)" "21.1% (259)" "19.9% (245)" "20.1% (247)"
paste(100 * round(table(dat1[dat1$LIT_grp=='yes',]$BL_ses)/dim(dat1[dat1$LIT_grp=='yes',])[1],3),
      "% (",
      table(dat1[dat1$LIT_grp=='yes',]$BL_ses), ")",sep='')

# control 
# ses
table(dat1[dat1$LIT_grp=='no',]$BL_ses)
100 * round(table(dat1[dat1$LIT_grp=='no',]$BL_ses)/dim(dat1[dat1$LIT_grp=='no',])[1],3)
#  "26.3% (325)" "21.1% (261)" "17.7% (219)" "18.4% (227)" "16.4% (203)"
paste(100 * round(table(dat1[dat1$LIT_grp=='no',]$BL_ses)/dim(dat1[dat1$LIT_grp=='no',])[1],3),
      "% (",
      table(dat1[dat1$LIT_grp=='no',]$BL_ses), ")",sep='')

col1 = c(col1,c('Poorest',
                'Poor',
                'Median poor',
                'Less poor',
                'Least poor'))
col2 = c(col2, paste(100 * round(table(dat1[dat1$LIT_grp=='yes',]$BL_ses)/dim(dat1[dat1$LIT_grp=='yes',])[1],3),
                     "% (",
                     table(dat1[dat1$LIT_grp=='yes',]$BL_ses), ")",sep=''))
col3 = c(col3, paste(100 * round(table(dat1[dat1$LIT_grp=='no',]$BL_ses)/dim(dat1[dat1$LIT_grp=='no',])[1],3),
                     "% (",
                     table(dat1[dat1$LIT_grp=='no',]$BL_ses), ")",sep=''))

# intervnetion
# baseline literacy – spelling score (0-20) – mean (sd)
round(mean(dat1[dat1$LIT_grp=='yes',]$BL_gll21_total),1)
round(sd(dat1[dat1$LIT_grp=='yes',]$BL_gll21_total),1)
# cell 8.4 (4.6)
paste(round(mean(dat1[dat1$LIT_grp=='yes',]$BL_gll21_total),1),
      ' (',
      round(sd(dat1[dat1$LIT_grp=='yes',]$BL_gll21_total),1),
      ')',sep='')

# control
# baseline literacy – spelling score (0-20) – mean (sd)
round(mean(dat1[dat1$LIT_grp=='no',]$BL_gll21_total),1)
round(sd(dat1[dat1$LIT_grp=='no',]$BL_gll21_total),1)
# 7.8 (4.3)
paste(round(mean(dat1[dat1$LIT_grp=='no',]$BL_gll21_total),1),
      ' (',
      round(sd(dat1[dat1$LIT_grp=='no',]$BL_gll21_total),1),
      ')',sep='')

col1 = c(col1,'Baseline literacy – spelling score (0-20) – mean (sd)')
col2 = c(col2, paste(round(mean(dat1[dat1$LIT_grp=='yes',]$BL_gll21_total),1),
                     ' (',
                     round(sd(dat1[dat1$LIT_grp=='yes',]$BL_gll21_total),1),
                     ')',sep=''))
col3 = c(col3, paste(round(mean(dat1[dat1$LIT_grp=='no',]$BL_gll21_total),1),
                     ' (',
                     round(sd(dat1[dat1$LIT_grp=='no',]$BL_gll21_total),1),
                     ')',sep=''))

### Outcome at 9-month follow-up 

col1 = c(col1, 'Outcome at 9-month follow-up ')
col2 = c(col2,' ')
col3 = c(col3,' ')

# intervention, outcome at 9-month follow-up 
sum(dat1[dat1$LIT_grp=='yes',]$gll21_total > 10,na.rm = TRUE)
100 * round(sum(dat1[dat1$LIT_grp=='yes',]$gll21_total > 10,na.rm = TRUE)/
              dim(dat1[dat1$LIT_grp=='yes', ])[1],3)
# cell "52.4% (644)"
paste(100 * round(sum(dat1[dat1$LIT_grp=='yes',]$gll21_total > 10,na.rm = TRUE)/
                    dim(dat1[dat1$LIT_grp=='yes', ])[1],3),
      "% (",
      sum(dat1[dat1$LIT_grp=='yes',]$gll21_total > 10,na.rm = TRUE),")",sep='')

# control
sum(dat1[dat1$LIT_grp=='no',]$gll21_total > 10,na.rm = TRUE)
100 * round(sum(dat1[dat1$LIT_grp=='no',]$gll21_total > 10,na.rm = TRUE)/
              dim(dat1[dat1$LIT_grp=='no', ])[1],3)
# cell "39.5% (488)"
paste(100 * round(sum(dat1[dat1$LIT_grp=='no',]$gll21_total > 10,na.rm = TRUE)/
                    dim(dat1[dat1$LIT_grp=='no', ])[1],3),
      "% (",
      sum(dat1[dat1$LIT_grp=='no',]$gll21_total > 10,na.rm = TRUE),")",sep='')

col1 = c(col1,'High literacy (spelling score > 10)')
col2 = c(col2, paste(100 * round(sum(dat1[dat1$LIT_grp=='yes',]$gll21_total > 10,na.rm = TRUE)/
                                   dim(dat1[dat1$LIT_grp=='yes', ])[1],3),
                     "% (",
                     sum(dat1[dat1$LIT_grp=='yes',]$gll21_total > 10,na.rm = TRUE),")",sep=''))
col3 = c(col3, paste(100 * round(sum(dat1[dat1$LIT_grp=='no',]$gll21_total > 10,na.rm = TRUE)/
                                   dim(dat1[dat1$LIT_grp=='no', ])[1],3),
                     "% (",
                     sum(dat1[dat1$LIT_grp=='no',]$gll21_total > 10,na.rm = TRUE),")",sep=''))

# intervention
# missing outcome
sum(is.na(dat1[dat1$LIT_grp=='yes',]$gll21_total))
100 * round(sum(is.na(dat1[dat1$LIT_grp=='yes',]$gll21_total))/ 
              dim(dat1[dat1$LIT_grp=='yes',])[1],3)
# cell 12.4% (152)
paste(100 * round(sum(is.na(dat1[dat1$LIT_grp=='yes',]$gll21_total))/ 
                    dim(dat1[dat1$LIT_grp=='yes',])[1],3),
      '% (',
      sum(is.na(dat1[dat1$LIT_grp=='yes',]$gll21_total)), ")",sep='')

# control
# missing outcome
sum(is.na(dat1[dat1$LIT_grp=='no',]$gll21_total))
100 * round(sum(is.na(dat1[dat1$LIT_grp=='no',]$gll21_total))/ 
              dim(dat1[dat1$LIT_grp=='no',])[1],3)
# cell 11.6% (143)
paste(100 * round(sum(is.na(dat1[dat1$LIT_grp=='no',]$gll21_total))/ 
                    dim(dat1[dat1$LIT_grp=='no',])[1],3),
      '% (',
      sum(is.na(dat1[dat1$LIT_grp=='no',]$gll21_total)), ")",sep='')
col1 = c(col1,'Missing outcome')
col2 = c(col2, paste(100 * round(sum(is.na(dat1[dat1$LIT_grp=='yes',]$gll21_total))/ 
                                   dim(dat1[dat1$LIT_grp=='yes',])[1],3),
                     '% (',
                     sum(is.na(dat1[dat1$LIT_grp=='yes',]$gll21_total)), ")",sep=''))
col3 = c(col3, paste(100 * round(sum(is.na(dat1[dat1$LIT_grp=='no',]$gll21_total))/ 
                                   dim(dat1[dat1$LIT_grp=='no',])[1],3),
                     '% (',
                     sum(is.na(dat1[dat1$LIT_grp=='no',]$gll21_total)), ")",sep=''))

### Combine the results to draw table:
table_1a = data.frame(col1, col2, col3)
colnames(table_1a) = c('',
                       'Intervention\n(n=1230)',
                       'Control\n(n=1235)')


table_1a
write.csv(table_1a,'table_1a.csv')




