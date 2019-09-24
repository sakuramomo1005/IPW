$x$

$ x^2 $
  


\begin{equation}
x^2
\end{equation}

$ \sum_{\forall i}{x_i^{2}} $

### Folders

* The **Code** folder contains all the R codes that used in our simulation

* The **Date** folder contains all the results from the simulations, which are saved as RData file and can be loaded in R environment directly

* The **Files** folder contains the summary files of the results

* The file: **Example codes.ipynb** is a notebook of a codes demo to show how we generate one dataset and how we analyze it




### Outcome model 

The outcome $Y_{ijl}$ is generated by:
$$\pi_{ijl}=expit(1+ 1.36 * i +  x_{ijl}+\delta_{ij}) $$ 

1. $i$ is the treatment indicator. $i=1$ treated; $i=0$ control

2. $x$ is the covariate. $x \sim N(0,0.2)$ (The variance is 0.2. I choose a relatively small variance since I want to avoid non-convergence. However, this generated small differences between CRA and adjusted CRA. )

3. The variance of $\delta_{ij}$ changed based on different ICC (the ICC for the datasets):

* ICC=0.05, $\delta_{ij} \sim N(0, 0.173)$ (the variance is 0.173)

* ICC=0.1, $\delta_{ij} \sim N(0, 0.366)$ (the variance is 0.366)

* ICC=0.2, $\delta_{ij} \sim N(0, 0.823)$ (the variance is 0.823)
  
  
4. Number of clusters: **10, 25, 50 ** clusters in each intervention arm

5. Cluster size: the cluster size is from a Poisson distribution, which is not fixed. cluster size $\sim POI(50)$ 


### Missingness generation model

$$ logit(R_{ijl}=1|Y_{ij},X_{ij})= intercept +  i  +  X_{ijl} + \sigma_{ij} $$
1. The intercept is varied to make sure the misisng percentage is around 30%

2. $i$ is the treatment indicator. $i$=1 treated; $i$=0 control

3. $x$ is the covariate. $x \sim N(0,0.2)$ (The variance is 0.2)

4. The variance of $\sigma_{ij}$ represents the cluster effects in missingness. It changed based on different missingness ICC:

* $\sigma_{ij}=0$, the missing ICC=0, there is no cluster effects

* ICC=0.1, $\sigma_{ij} \sim N(0, 0.366)$ (the variance is 0.366)

* ICC=0.3, $\sigma_{ij} \sim N(0, 1.410)$ (the variance is 1.410)

* ICC=0.5, $\sigma_{ij} \sim N(0, 3.291)$ (the variance is 3.291)

* ICC=0.7, $\sigma_{ij} \sim N(0, 7.678)$ (the variance is 7.678)

* ICC=0.9, $\sigma_{ij} \sim N(0, 29.616)$ (the variance is 29.616)

5. 1000 replicates for each scenario         

### Missingness handling methods

**1. Calculation of true value**: 

* Fit the GEE with the formula: Y $\sim$ intervention arm (without covariates). Estimate the coefficient of intervention arm.  

* With full datasets without missing values

* Repeat for 1000 times and calculate the mean value. 


**2. UCRA: unadjusted complete record analysis**

* Fit the GEE with the formula: Y $\sim$ intervention arm (without covariates). Estimate the coefficient of intervention arm.  

* Delete the records with missing values in Y

* Repeat for 1000 times. 


**3. CRA: adjusted complete record analysis**
 
* Fit the GEE with formula: Y $\sim$ intervention-arm + X (with covariates). Estimate the coefficient of intervention arm.  

* Delete the records with missing values in Y

* Repeat for 1000 times. 


**4. IPW: inverse probability weighting**

* Calculate the weights by fitting GLM: glm(y $\sim$ arm + x)

* Fit the GEE with the formula: Y $\sim$ intervention arm (without covariates) with corresponding weights. Estimate the coefficient of intervention arm.  

* Repeat for 1000 times. 



**5. IPWC: inverse probability weighting with cluster effects**

* Calculate the weights by fitting generalized linear mixed effect model: glmer (y $\sim$ arm + x + cluster effect)

* Fit the GEE with the formula: Y $\sim$ intervention arm (without covariates) with corresponding weights. Estimate the coefficient of intervention arm.  

* Repeat for 1000 times. 

**6. MMI: multilevel multiple imputation**

* consider cluster effects in the imputation

