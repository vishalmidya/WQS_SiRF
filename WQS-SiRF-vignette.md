# How to implement the WQS-SiRF algorithm
## Vishal Midya and Chris Gennings

This article presents a step-by-step guide and intuitions to implement the WQS-SiRF algorithm. WQS-SiRF is conducted in two stages to discover interactions associated with the outcome of interest. The first part of this algorithm uses an exposure mixture model (which, in this case, is the Weighted Quantile Sum regression (WQS) model). The WQS mixture model was used because of the apriori hypothesized directionality of the mixture association. Note that other exposure mixtures analytical models, such as the Bayesian Kernel Machine Regression (BKMR) or Quantile g-computation method, can also be used when the directionality of the association between exposures and the outcome is not hypothesized beforehand, and interest lies in the overall mixture effect. See more technical details in Midya et al. <sup>1</sup>.

Assuming the main chemical mixture effect and the interactions are additive; we extracted the residuals (Pearson residuals in case of binary outcome) from this model. We treated the residual as the new outcome and used a machine learning-based prediction framework to discover non-linear combinations predictive of the outcome. The interactions are searched using a repeated hold-out signed-iterated Random Forest (rh-SiRF), where the predictors are concentrations of the chemical exposures. The SiRF (Signed Iterative Random Forest) algorithm combined a state-of-the-art predictive tool called "Iterative Random Forests" with "Random Intersection Trees" to search for combinations predictive of the outcome <sup>1</sup>. On top of the SIRF algorithm, we introduced a repeated hold-out step that randomly partitions the data in training and testing sets for better generalizability. See more technical details in Basu et al.<sup>2</sup>.

As noted in Midya et al.<sup>1</sup>, this combination of exposure mixture algorithms can search for multi-ordered and non-linear interactions even when the number of chemical exposures is large. Further, since these interactions are based on thresholds, i.e., these interactions mimic the classical toxicological paradigm in which an interaction occurs only if the concentrations of certain chemicals are above some threshold. 


##  Simulated exposure data

We first present simulated data on around `500` hypothetical participants with `60` simulated Taxa. The dataset, named `data.simulated.csv`, is uploaded as a part of the demonstration.

## Required `R` packages


Kindly install the following `R` packages
1. `mvtnorm`
2. `iRF`

Find the instructions [here](https://github.com/sumbose/iRF/tree/master) to install the `iRF` package.
