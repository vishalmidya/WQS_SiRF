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

## Data generating Process: creating a 3<sup>rd</sup> ordered non-linear interaction and an outcome

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
All the `25` exposures are already supplied in the dataset. These are named, `V1`, `V2`,..., and `V25`.  

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
This outcome is already supplied in the dataset. However, note that the effect size of the `three_clique` is set at `1` and the mixture effect is set at `0.5`

## The aim of this algorithm

1. To recover the mixture effect (joint mixture effect in apriori selected direction - WQS or overall mixture effect - qgcomputation)
2. To recover the 3<sup>rd</sup> order interaction of `V1`, `V3`, and `V11`.
3. To recover the thresholds for the construction of the interaction.
4. To estimate the effect of the interaction on the outcome.

Note that, there are a total of `300` two-ordered, and `2300` three-ordered combinations to choose from `25` exposures. On top of that, one needs to find the thresholds and determine the directionality.  Further, as the number of exposures increases, the total number of possible combinations to test from exponentially increases. Therefore, we used the signed-iterated Random Forest (SiRF) algorithm to search for the optimal combinations. Further, we introduced a repeated-holdout stage on the SiRF algorithm to keep the false positives and negatives in check.

## WQS-SiRF Algorithm

1. __Generalized Weighted Quantile Sum Regression__

First, we fitted a `gwqs` model. We implemented the random subset and repeated holdout variants of WQS <sup>3,4</sup> to make the analysis robust. We kept `40%` of the data for validation, and repeated the holding-out procedure for `50` iterations, with `200` bootstraps each time. Given the apriori hypothesized directionality (think the adverse effect of chemical mixture), we chose `b1_pos = T`. 

```{}
exposures <- paste0("V",seq(1:25))
gg <- gwqsrh(outcome ~ wqs + cov1 + cov2 + cov3 + cov4,
           mix_name = exposures, data = data.simulated, q = NULL, signal = "t2",
           validation = 0.4, b = 200, rs = T,  n_vars = p0/2, rh = 50, 
           b1_pos = T,  family = "gaussian",plan_strategy = "multicore")

summary(gg)
```
The result is shown below:

```{}
             Estimate Std. Error   2.5 % 97.5 %
 (Intercept)  0.17661    0.01873 0.13990  0.213
 wqs          0.61232    0.04246 0.52910  0.696
 cov1         0.07110    0.02669 0.01878  0.123
 cov2         0.07065    0.02662 0.01847  0.123
 cov3         0.07263    0.02967 0.01448  0.131
 cov4         0.06899    0.02151 0.02684  0.111
```

2. __Extract the residuals from the Generalized Weighted Quantile Sum Regression__

The extracted residual from this `gwqs` model, named `wqs.residuals`, has been added to the simulated dataset. Note that, for a binary outcome, one can use the _Pearson residuals_, which exhibit asymptotic normal properties for large samples. This residual, `wqs.residuals` is the outcome.

Further, note that individually `(wqs.residual ~ exposure + cov1 + cov2 + cov3 + cov4)`, the concentrations of the three exposures may _not_ be significantly associated with the outcome. But through the three-ordered interaction, there is a significant statistical association.  

```{}
             Estimate Std. Error t value  Pr(>|t|)
V1          -0.061736   0.020194  -3.057  0.00235 **
V3           0.005712   0.019172   0.298    0.766
V11          0.025359   0.020004   1.268    0.205
three_clique  0.87192    0.03281   26.57   <2e-16 ***

```







References:

1. ff
2. ff
3. Tanner, E. M.; Bornehag, C.-G.; Gennings, C. Repeated holdout validation for weighted quantile sum regression. MethodsX. 2019, 6, 2855−2860.
4. Curtin, P.; Kellogg, J.; Cech, N.; Gennings, C. A random subset implementation of weighted quantile sum (WQSRS) regression for analysis of high-dimensional mixtures. Communications in Statistics - Simulation and Computation. 2021, 50, 1119−1134.


