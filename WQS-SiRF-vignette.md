# How to implement the WQS-SiRF algorithm
## Vishal Midya and Chris Gennings

This article presents a step-by-step guide and intuitions to implement the WQS-SiRF algorithm. WQS-SiRF is conducted in two stages to discover interactions associated with the outcome of interest. The first part of this algorithm uses an exposure mixture model (which, in this case, is the Weighted Quantile Sum regression (WQS) model). The WQS mixture model was used because of the apriori hypothesized directionality of the mixture association. Note that other exposure mixtures analytical models, such as the Bayesian Kernel Machine Regression (BKMR) or Quantile g-computation method, can also be used when the directionality of the association between exposures and the outcome is not hypothesized beforehand, and interest lies in the overall mixture effect. See more technical details in Midya et al. <sup>1</sup>.

Assuming the main chemical mixture effect and the interactions are additive; we extracted the residuals (Pearson residuals in case of binary outcome) from this model. We treated the residual as the new outcome and used a machine learning-based prediction framework to discover non-linear combinations predictive of the outcome. The interactions are searched using a repeated hold-out signed-iterated Random Forest (rh-SiRF), where the predictors are concentrations of the chemical exposures. The SiRF (Signed Iterative Random Forest) algorithm combined a state-of-the-art predictive tool called "Iterative Random Forests" with "Random Intersection Trees" to search for combinations predictive of the outcome <sup>1</sup>. On top of the SIRF algorithm, we introduced a repeated hold-out step that randomly partitions the data in training and testing sets for better generalizability. See more technical details in Basu et al.<sup>2</sup>.

As noted in Midya et al.<sup>1</sup>, this combination of exposure mixture algorithms can search for multi-ordered and non-linear interactions even when the number of chemical exposures is large. Further, since these interactions are based on thresholds, i.e., these interactions mimic the classical toxicological paradigm in which an interaction occurs only if the concentrations of certain chemicals are above some threshold. 

The following illustration included both the _Weighted Quantile Sum regression (WQS) regression_ and _Quantile g-computation_ models in conjunction with the rh-SiRF algorithm to find a three-ordered non-linear interaction. 

##  Simulated exposure data

We first present simulated data on around `500` hypothetical participants with `60` simulated V The dataset, named `data.simulated.csv`, is uploaded as a part of the demonstration.

## Required `R` packages


Kindly install the following `R` packages
1. `mvtnorm`
2. `gWQS`
3. `iRF`

Find the instructions [here](https://github.com/sumbose/iRF/tree/master) to install the `iRF` package.

See details on how to use the `gWQS` model [here](https://cran.r-project.org/web/packages/gWQS/vignettes/gwqs-vignette.html).

## Creating a 3<sup>rd</sup> ordered non-linear interaction and an outcome

1. __Run the following chunk of functions:__

`require(mvtnorm)`

```{}
make.SIGMA <- function(rho,dimension){
  SIGMA = matrix(NA,dimension,dimension)
  for(i in 1:dimension){
    for(j in 1:dimension){
      if(i != j){
        a <- sign(rnorm(1,0,1))*rnorm(1,0.1, 0.01)
        if(a>0) {SIGMA[i,j] = rho + a}
        else if (a <0)  {SIGMA[i,j] = rho}
      }
      else if (i == j){
        SIGMA[i,j] = 1
      }
    }
  }
  
  (SIGMA + t(SIGMA))/2
}
```

```{}
make.X0 <- function(n, rho, p0){
  X0 <- (mvtnorm::rmvnorm(n, mean = rep(0, nrow(make.SIGMA(rho, p0))), sigma = make.SIGMA(rho, p0)))
}
```
The set of `25` correlated exposures was created using the following code.

```{}
n = 500
p0 = 25
set.seed(123456)
X <- make.X0(n, 0.2, p0)
```
All the `25` exposures are already supplied in the dataset. These are named as, `V1`, `V2`,..., and `V25`.  

2. __Next, we create four covariates that are moderately correlated to each other__

```{}
n <- dim(sample.data.simulated)[1]
set.seed(12456)
covariates <- make.X0(n,0.1,4)
```
The set of four covariates is already supplied in the dataset.  

3. __Create a 3<sup>rd</sup> ordered non-linear interaction__

Without loss of generality, we choose 
  - `V1`, `V5`, `V10`, `V15`, and `V25` to form a positive mixture effect associated with the outcome.
  
  - On top of the mixture effect, `V1`, `V3`, and `V11` to form a three-ordered interaction. This interaction is present only in those individuals where (1) the concentration of `V1` is _less_ than its 60<sup>th</sup> percentile, (2) the concentration of `V3` is _greater_ than its 40<sup>th</sup> percentile, and lastly, (3) the concentration of `V11` is _greater_ than its 40<sup>th</sup> percentile. In terms of code,

```{}
three_clique <- as.numeric(data.simulated$V1 <= quantile(data.simulated$V1, 0.6)
                             & data.simulated$V3 >= quantile(data.simulated$V3, 0.4)
                                & data.simulated$V11 >= quantile(data.simulated$V11, 0.4))
```

4. __Create the outcome__

We create a Gaussian outcome composed of the mixture effect and the indicator for the three-ordered interaction, covariates, and random Gaussian error.

```{}
set.seed(45667)
(X[,c(1,5,10,15,20)] %*% c(0.1,0.1, 0.1, 0.1,0.1)) +  three_clique + covariates %*% c(0.06,0.06,0.06,0.06) + rnorm(n, 0, 0.2)
```
This outcome is already supplied in the dataset. However, note that the effect size of the `three_clique` is set at `1`. Further, note that individually `(outcome ~ Taxa + cov1 + cov2 + cov3 + cov4)`, the relative abundances of the three Taxa are _not_ significantly associated with the outcome. Only through the three-ordered clique, is there a significant statistical association.  

```{}
             Estimate Std. Error t value Pr(>|t|)
V1      -0.009965   0.006567  -1.517 0.129816    
V3       0.007614   0.006652   1.145  0.25297    
V11      0.008184   0.006479   1.263 0.207165    
clique        1.24607    0.06482  19.223  < 2e-16 ***

```

