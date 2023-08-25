# How to implement the WQS-SiRF algorithm
## Vishal Midya and Chris Gennings

This article presents a step-by-step guide and intuitions to implement the WQS-SiRF algorithm. WQS-SiRF is conducted in two stages to discover interactions associated with the outcome of interest. The first part of this algorithm uses an exposure mixture model (which, in this case, is the Weighted Quantile Sum regression (WQS) model). The WQS mixture model was used because of the apriori hypothesized directionality of the mixture association. Note that other exposure mixtures analytical models, such as the Bayesian Kernel Machine Regression (BKMR) or Quantile g-computation method, can also be used when the directionality of the association between exposures and the outcome is not hypothesized beforehand, and interest lies in the overall mixture effect. See more technical details in Midya et al. <sup>1</sup>.

Assuming the main chemical mixture effect and the interactions are additive; we extracted the residuals (Pearson residuals in case of binary outcome) from this model. We treated the residual as the new outcome and used a machine learning-based prediction framework to discover non-linear combinations predictive of the outcome. The interactions are searched using a repeated hold-out signed-iterated Random Forest (rh-SiRF), where the predictors are concentrations of the chemical exposures. The SiRF (Signed Iterative Random Forest) algorithm combined a state-of-the-art predictive tool called "Iterative Random Forests" with "Random Intersection Trees" to search for combinations predictive of the outcome <sup>1</sup>. On top of the SIRF algorithm, we introduced a repeated hold-out step that randomly partitions the data in training and testing sets for better generalizability. See more technical details in Basu et al.<sup>2</sup>.

As noted in Midya et al.<sup>1</sup>, this combination of exposure mixture algorithms can search for multi-ordered and non-linear interactions even when the number of chemical exposures is large. Further, since these interactions are based on thresholds, i.e., these interactions mimic the classical toxicological paradigm in which an interaction occurs only if the concentrations of certain chemicals are above some threshold. 

The following illustration included the _Weighted Quantile Sum regression (WQS) regression_ in conjunction with the rh-SiRF algorithm to find a three-ordered non-linear interaction. However, one can replace the WQS model with the _Quantile g-computation_ model, but caution should be practiced while interpreting the result (more details later). 

##  Simulated exposure data

We first present simulated data on around `500` hypothetical participants with `25` simulated exposures. The dataset, named `data.simulated.csv`, is uploaded as a part of the demonstration.

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

The extracted residual from this `gwqs` model, named `wqs.residuals`, has been added to the simulated dataset. Note that, for a binary outcome, one can use the _Pearson residuals_, which exhibit asymptotic normal properties for large samples. This residual, `wqs.residuals` is the outcome. Further, note that individually `(wqs.residual ~ exposure + cov1 + cov2 + cov3 + cov4)`, the concentrations of the three exposures may _not_ be significantly associated with the outcome. But through the three-ordered interaction, there is a significant statistical association.  

```{}
             Estimate Std. Error t value  Pr(>|t|)
V1          -0.061736   0.020194  -3.057  0.00235 **
V3           0.005712   0.019172   0.298    0.766
V11          0.025359   0.020004   1.268    0.205
three_clique  0.87192    0.03281   26.57   <2e-16 ***

```
Note that, there are a total of `300` two-ordered, and `2300` three-ordered combinations to choose from `25` exposures. On top of that, one needs to find the thresholds and determine the directionality.  Further, as the number of exposures increases, the total number of possible combinations to test from exponentially increases. Therefore, we used the signed-iterated Random Forest (SiRF) algorithm to search for the optimal combinations. Further, we introduced a repeated-holdout stage on the SiRF algorithm to keep the false positives and negatives in check.

## Finding the optimal combination of exposures

Run the following function that finds the most frequently occurring combinations of exposures:
Copy this code chunk below and run it.

```{}
require("iRF")
clique.finder <- function(exposures, outcome, iterations, validation, seed.value, n.bootstrap, min.prevalence, min.stability, data){
  n <- dim(data)[1]
  if(n == 0){
    stop("Provide a dataset in data.frame format") 
  }
  if(length(exposures) < 2){
    stop("Need at least two exposures to find a clique") 
  }
  if(validation == 0 | validation >=1){
    stop("Validation must be a positive fraction within 0 and 1") 
  }
  train.sample <- 1 - validation
  if(iterations < 50){
    warning("Use more than 50 iterations for reliable replication")
  }
  if(n.bootstrap < 100){
    warning("Use more than 100 bootstraps for reliable replication")
  }
  tab <- data.frame( int = NA_character_,  prevalence = NA_real_, precision = NA_real_, cpe = NA_real_,
                     sta.cpe = NA_real_,  fsd = NA_real_, sta.fsd = NA_real_, mip = NA_real_,
                     sta.mip = NA_real_, stability =NA_real_, occurance.freq = NA_character_ , occurance.freq.sum = NA_real_)
  for(i in 1:iterations){
    set.seed(runif(1, 0, (seed.value + 10)))
    train.id <- sample(seq(1,n), ceiling(n*train.sample))
    test.id <- setdiff(1:n, train.id)
    fit <- iRF(x=data.simulated[train.id, exposures], 
               y=data.simulated[train.id,outcome], 
               xtest=data.simulated[test.id, exposures],
               ytest=data.simulated[test.id,outcome],
               n.iter=10, 
               n.core=3,
               select.iter = T,
               n.bootstrap=n.bootstrap
    )
    SiRF_table <- as.data.frame(fit$interaction[fit$interaction$stability > min.stability,])
    if(dim(SiRF_table)[1] != 0){
      
      SiRF_table_subset <- SiRF_table
      SiRF_table_subset_prev <- SiRF_table_subset[SiRF_table_subset$prevalence > min.prevalence,]
      
      if(dim(SiRF_table_subset_prev)[1] != 0){
        nam <- NA_character_
        for(i in 1:nrow(SiRF_table_subset_prev)){nam <- c(nam, strsplit(SiRF_table_subset_prev$int[i],"_")[[1]])}
        yt<-as.data.frame(table(nam))
        SiRF_table_subset_prev$occurance.freq <- rep(NA_character_,nrow(SiRF_table_subset_prev))
        for(i in 1:nrow(SiRF_table_subset_prev)){
          SiRF_table_subset_prev$occurance.freq[i] <- as.character((paste0(yt$Freq[which(yt$nam == strsplit(SiRF_table_subset_prev$int[i],"_")[[1]][1])],"/",yt$Freq[which(yt$nam == strsplit(SiRF_table_subset_prev$int[i],"_")[[1]][2])])))
        }
        SiRF_table_subset_prev$occurance.freq.sum <- rep(NA_character_,nrow(SiRF_table_subset_prev))
        for(i in 1:nrow(SiRF_table_subset_prev)){
          SiRF_table_subset_prev$occurance.freq.sum[i] <- sum(strsplit(SiRF_table_subset_prev$occurance.freq[i],"/")[[1]] %in% "1")
        }
        tab <- rbind(tab, SiRF_table_subset_prev)
      } else {
        tab <- rbind(tab, SiRF_table_subset_prev)
      }
    } else {
      tab <- rbind(tab, SiRF_table)
    }
  }
  clique <- tab
  if(dim(tab)[1] == 1){
    stop("No clique was found, decrease the value of min.stability and increase the number of bootstraps")
  }
  clique <- clique[-1,]
  for(i in 1:nrow(clique)){
    if(sum(strsplit(clique$int[i],"")[[1]] %in% "+") != 0 & sum(strsplit(clique$int[i],"")[[1]] %in% "-") == 0){
      x <- strsplit(clique$int[i],"+")[[1]]
      clique$int[i] <- paste0(x[x!= "+"], collapse = "")  
    } else if(sum(strsplit(clique$int[i],"")[[1]] %in% "-") != 0 & sum(strsplit(clique$int[i],"")[[1]] %in% "+") == 0){
      x <- strsplit(clique$int[i],"-")[[1]]
      clique$int[i] <- paste0(x[x!= "-"], collapse = "")  
    } else if (sum(strsplit(clique$int[i],"")[[1]] %in% "-") != 0 & sum(strsplit(clique$int[i],"")[[1]] %in% "+") != 0){
      x <- strsplit(clique$int[i],"")[[1]]
      clique$int[i] <- paste0(x[x!= "-" & x!= "+"], collapse = "")  
    }
  }
  rtf <- as.data.frame(table(clique$int))
  rtf <- rtf[order(rtf$Freq, decreasing = T),]
  return(as.data.frame(rtf))
}
```

Finally, run the `function` called `clique.finder`. This is the main function where you can change its arguments. Below we discuss each argument for this function and what they entail.

```{}
clique.finder(exposures = paste0("V", seq(1,25)), outcome = "wqs.residuals",  iterations = 500, validation = 0.4, 
              seed.value = 1234, n.bootstrap = 200, min.prevalence = 0.05, min.stability = 0.25, data = data.simulated)
```

1. `exposures`: a vector of all possible exposures (among which one intends to find the combinations)
2. `outcome`: name of the outcome variable
3. `iterations`: the number of repeated holdouts (should be more than 100)
4. `validation`: the proportion of the dataset which is set aside for validation at each repeated holdout iteration 
5. `seed.value`: random initial seed value for partitioning the dataset
6. `n.bootstrap`: the number of bootstrap iterations employed at each of the repeated holdouts on the training dataset (should be more than 100)
7. `min.prevalence`: the minimum proportion (lower bound) of the sample that has the clique. Here we chose `5%` as the lower bound of the prevalence. 
8. `min.stability`: the stability implies the proportion of times the combination was recovered across bootstrap replicates. The `min.stability` is the lower bound. Here we chose `25%` as the lower bound.
9. `data`: name of the dataset

Note that lowering the values of `min.prevalence` and `min.stability` finds more combinations of exposures; however, due to the implementation of the repeated holdout technique, lowering these bounds does not significantly affect the most stable or most frequently occurring combinations. In this example, we used `500` repeated holdouts. One can utilize the [parallel R package](https://www.rdocumentation.org/packages/parallel/versions/3.6.2) for fast parallel computation. The same code can be used for binary outcomes as well.

Here is the result of the top 10 combinations obtained from the simulated dataset

```{}
         Var1    Freq
       V11_V3 19.1023
        V1_V3 18.9960
       V1_V11 16.8688
    V1_V11_V3  5.3499
       V1_V19  3.8609
       V1_V25  3.2972
       V1_V24  3.0951
        V1_V8  2.4782
       V25_V3  2.3931
        V1_V7  2.3825
```
The first column, `Var1`, denotes all possible combinations picked up by the algorithm, and the `Freq` denotes the percentage of their occurrence over all possible detected combinations and all bootstrap and repeated holdout combinations. The total number of detected unique combinations is around `80`. Note the top three combinations, `V11/V3`, `V1/V11`, and `V1/V3`, occurred the most (more than 10%) among all possible detected combinations. Moreover, a three-ordered `V1/V11/V3` combination was among the most occurring. The three-ordered combination of `V1/V11/V3` induced a downstream of further two-ordered combinations. With smaller sample sizes, this tool effectively detects lower-ordered combinations. However, as the sample size increases, the frequency of the true higher-ordered combinations also increases. Although this is very subjective, a good practice is to choose the top (first three or first five) most frequently occurring combinations as long as they form a clique.  Although we found the exposure combination, how it is associated with the outcome or the directionality is unknown. In the next stage, we estimate the thresholds of the exposures and their joint association with the outcome. 

## Estimating the thresholds for the interacting exposures and the joint association with the outcome

Run the following function that finds the thresholds for the exposures `V1`, `V3`, and `V11`. Each exposure's directionality is chosen based on its univariate association with the outcome. The following code _should only be used_ based on the output from the `clique.finder` function; otherwise, overfitting is possible. 

```{}

clique.tba <- function(clique.names, outcome, grid.quantile, min.prevalence,  data, family){
  len <- length(clique.names)
  if(len < 2){
    stop("Need at least two exposures to form a meaningful clique") 
  }
  n <- dim(data)[1]
  if(n == 0){
    stop("Please provide a dataset in data.frame format") 
  }
  beta.data <- data.frame(Exposure = rep(NA_character_, len), effect_size = rep(NA_real_, len))
  beta.data$Exposure <- clique.names
  for(i in 1:len){
    g.out <- data[,outcome]
    if(family == "gaussian"){
      fit <- summary(lm(g.out ~ as.matrix(data[, c(beta.data$Exposure[i])]), data = data))  
    }
    if(family == "binomial"){
      fit <- summary(glm(g.out ~ as.matrix(data[, c(beta.data$Exposure[i])]), data = data , family = family))  
    }
    if(family == "poisson"){
      fit <- summary(glm(g.out ~ as.matrix(data[, c(beta.data$Exposure[i])]), data = data , family = family))  
    }
    beta.data$effect_size[i] <- fit$coefficients[2,1]
  }
  x <- grid.quantile
  if(length(grid.quantile) == 1){
    stop("Please increase the number of possible thresholds")
  }
  if(sum(grid.quantile >= 1) != 0){
    stop("All the values provided must be less than 1")
  }
  if(sum(grid.quantile <= 0) != 0){
    stop("All the values provided must be more than 0")
  }
  d1 <- do.call(expand.grid, replicate(len, x, simplify = F))
  d1$min.prevalence <- rep(NA_real_, dim(d1)[1])
  d1$effect_size <- rep(NA_real_, dim(d1)[1])
  d1$se <- rep(NA_real_, dim(d1)[1])
  d1$pvalue <- rep(NA_real_, dim(d1)[1])
  for(i in 1:nrow(d1)){
    mat.len <- as.data.frame(matrix(NA_real_, ncol = len, nrow = nrow(data)))
    for(j in 1:len){
      if(sign(beta.data$effect_size)[j] < 0){
        mat.len[,j] <- as.numeric(data[,clique.names[j]] <= quantile(data[,clique.names[j]], d1[i,j] )) 
      }else{
        mat.len[,j] <- as.numeric(data[,clique.names[j]] >= quantile(data[,clique.names[j]], d1[i,j] ))
      }
    }
    clique.int <- apply(mat.len, 1, function(x){prod(x)})
    
    if(sum(clique.int)!= 0){
      
      d1$min.prevalence[i] <- as.numeric(table(clique.int)/sum(table(clique.int)))[2]
      
      if(d1$min.prevalence[i] >= min.prevalence){
        
        data$clique.int <- clique.int
        g.out <- data[,outcome]
        if(family == "poisson"){
          s <- summary(glm(g.out ~ as.matrix(data[, c("clique.int")]), data = data, family = family))
        }
        if(family == "binomial"){
          s <- summary(glm(g.out ~ as.matrix(data[, c("clique.int")]), data = data, family = family))
        }
        if(family == "gaussian"){
          s <- summary(lm(g.out ~ as.matrix(data[, c("clique.int")]), data = data))
        }
        d1$effect_size[i] = s$coefficients[2,1]
        d1$se[i] = s$coefficients[2,2]
        d1$pvalue[i] = s$coefficients[2,4]
      }
    }
  }
  d2 <- na.omit(d1)
  if(dim(d2)[1] == 0){
    stop("Please decrease the min.prevalence, but keep it more than 5% for reliability")
  }
  d2 <- d2[d2$min.prevalence > 0.1,]
  d2 <- d2[order(abs(d2$effect_size), decreasing = T),]
  out <- d2[1,]
  
  cuts <-colnames(out)[!(colnames(out) %in% c("min.prevalence", "effect_size", "se",  "pvalue" ))]
  colnames(out)[1:length(cuts)] <- paste0(clique.names, ":Threshold")
  for(i in 1:nrow(beta.data)){
    if(beta.data[i,"effect_size"] < 0){
      out[,i] <- paste0("<=", out[,i]*100,"th Percentile")
    } else {out[,i] <- paste0(">=", out[,i]*100,"th Percentile")}
  }
  return(out)
}

```
Finally, run the `function` called `clique.tba`. Below we discuss each argument for this function and what they entail.

```{}
clique.tba(clique.names = c("V1", "V3", "V11"), outcome= "wqs.residuals", 
           grid.quantile = seq(0.2, 0.8, 0.1), min.prevalence = 0.05, family = "gaussian", data = data.simulated)
```
1. `clique.names`: a vector containing the names of the most frequently occurring exposures obtained from the previous code chunk. We chose `V1`, `V3`, and `V11` based on the output from `clique.finder`.
2. `outcome`: name of the outcome variable
3. `grid.quantile`: choices of the quantiles to search for the thresholds. We removed the lower and upper 20<sup>th</sup> quantiles for stable results. 
4. `min.prevalence`: the minimum proportion (lower bound) of the sample that has the combination. Here we chose `10%` as the lower bound of the prevalence. 
5. `family`: choice of `glm` family of distributions
6. `data`: name of the dataset

I'm sharing below the final output from the simulated example.

```{}
      V1:Threshold      V3:Threshold     V11:Threshold min.prevalence effect_size         se       pvalue
 <=60th Percentile >=40th Percentile >=40th Percentile          0.188   0.8719229 0.03281273 1.558698e-97
```

The `V1:Threshold`, `V3:Threshold`, and `V11:Threshold` denote the estimated thresholds for `V1`, `V3`, and `V11`, respectively. Therefore, this exposure combination is formed in those having (1) `V1` less than 60<sup>th</sup> percentile of the sample, (2) `V3` greater than 40<sup>th</sup> percentile of the sample, and lastly (3) `V11` more than 40<sup>th</sup> percentile of the sample. The recovered estimated effect size of this three-ordered interaction is `0.9`, and the estimated prevalence of this exposure combination is almost `19%`.  Note that, a similar result is obtained if one chooses to use the _Quantile g-computation_ on this simulated dataset, instead of the _Weighted Quantile Sum regression (WQS) regression_. This is the case since, in the simulation, all five of the exposures were positively associated with the outcome. Therefore, in the absence of any other effect, the overall mixture effect from Quantile g-computation is similar to the joint mixture effect in the positive direction.

### References:

1. Midya, V., Alcala, C. S., Rechtman, E., Gregory, J. K., Kannan, K., Hertz-Picciotto, I., ... & Valvi, D. (2023). Machine Learning Assisted Discovery of Interactions between Pesticides, Phthalates, Phenols, and Trace Elements in Child Neurodevelopment. Environmental Science & Technology.
2. Basu, S.; Kumbier, K.; Brown, J. B.; Yu, B. Iterative Random Forests to Discover Predictive and Stable High-Order Interactions. Proc. Natl. Acad. Sci. 2018, 115 (8), 1943–1948.
3. Tanner, E. M.; Bornehag, C.-G.; Gennings, C. Repeated holdout validation for weighted quantile sum regression. MethodsX. 2019, 6, 2855−2860.
4. Curtin, P.; Kellogg, J.; Cech, N.; Gennings, C. A random subset implementation of weighted quantile sum (WQSRS) regression for analysis of high-dimensional mixtures. Communications in Statistics - Simulation and Computation. 2021, 50, 1119−1134.

### Acknowledgments

This method was developed at the Dept. of Environmental Medicine and Public Health, Icahn School of Medicine at Mount Sinai, NYC, with funding and support from the National Institute of Environmental Health Sciences (P30ES023515, and U2C ES026555-01).

### Contact

Please email Vishal Midya (vishal.midya@mssm.edu) for reporting typos, general questions, and feedback.  
