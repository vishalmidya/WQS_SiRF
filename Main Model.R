library("xgboost")
library('gWQS')
library('MASS')
library("car")
library("mvtnorm")
library("dplyr")
library('doMC')
library('doParallel')
library("gbm")
library("iRF")
library('mice')
library('ggplot2')

# Import Targeted data

targetted_data <- read.csv("....\\1461_TARGETED_DATA.csv")
targetted_data$analyte_code[targetted_data$analyte_code == "TCP24" & targetted_data$LOD == 0.1] <- "TCP245"
targetted_data$analyte_code[targetted_data$analyte_code == "TCP24" & targetted_data$LOD == 0.02] <- "TCP246"
targetted_data$CHILD_PID <- as.character(targetted_data$CHILD_PID)

# Getting rid of PFAS analytes

targetted_data <- targetted_data[!(targetted_data$analyte_code %in% c('NETFOSAA', 'NMFOSAA', 'PFBS', 'PFDA', 'PFDODA', 'PFHPA', 'PFHXA', 'PFHXS', 'PFNA', 'PFOA', 'PFOS',  'PFOSA', 'PFPEA', 'PFUNDA')),]

# Adjusting for Specific Gravity

targetted_data$cf <- (1.0223 - 1)/(targetted_data$SG_avg - 1)
for(i in 1:nrow(targetted_data)){
  if(targetted_data$cf[i] <= 0.5){targetted_data$cf[i] = 0.5}
  else if (targetted_data$cf[i] >= 2){targetted_data$cf[i] = 2}
}

# Imputing the below LOD values by LOD/2

for(i in 1:nrow(targetted_data)){
  if(targetted_data$comment_code[i] == 37){
    targetted_data$concentration[i] = targetted_data$LOD[i]/2
  }
}

targetted_data$sgc_concentration <- targetted_data$cf*targetted_data$concentration

# Log (base = 2) transformation of concentrations

targetted_data$sgc_log_concentration <- (log(targetted_data$sgc_concentration, 2))

# Histogram of the specific gravity and log transformed concentration  

hist(targetted_data$sgc_log_concentration)


# Converting the data from Long to Wide format

concen <- as.data.frame(reshape(targetted_data[,c("CHILD_PID","analyte_code","sgc_log_concentration","sample_collection_year")], 
                                idvar = c("CHILD_PID","sample_collection_year"), 
                                timevar = "analyte_code", direction = "wide"))

df <- unlist(strsplit(colnames(concen)[c(-1,-2)],"sgc_log_concentration."))
colnames(concen) <- c("CHILD_PID","sample_collection_year",df[df != ""])


# Importing the Epi data

epi <- read.csv("....\\1461_EPI_DATA.csv")
d <- merge(epi, concen, by = "CHILD_PID")
df <- d[,colnames(d) %in% c(colnames(concen),"correctsex","YOB","ChildAgeMos","ChildRace_3cat","OwnHome","MMC_4cat","AgeMomYrs","dx2_group")]
colnames(df)[3] <- "csex"
colnames(df)[4] <- "crace"
colnames(df)[5] <- "cage"
df$cage <- df$cage/12
df$YOB <- df$YOB - 2000

# Imputing some missing covariates using MICE

init = mice(df[!(colnames(df) %in% c("CHILD_PID","dx2_group"))], maxit=0)
meth = init$method
predM = init$predictorMatrix
set.seed(1234)
imputed = mice(df[!(colnames(df) %in% c("CHILD_PID","dx2_group"))], method="pmm", predictorMatrix=predM, m=3, ntree = 100)
imputed <- complete(imputed, action = 3) # this is your final imputed dataset

finaldf <- df
finaldf$crace <- as.character(imputed$crace)
finaldf$MMC_4cat <- as.character(imputed$MMC_4cat)
finaldf$OwnHome <- as.character(imputed$OwnHome)

finaldf$ASD <- df$dx2_group
finaldf <- finaldf[finaldf$ASD ==2 | finaldf$ASD == 5,]
finaldf$ASD <- ifelse(finaldf$ASD == 2, 1, 0)

# The final dataset with exposures and the outcome is 'finaldf'

################################################################################

# Analysis

## Weighted Quantile Sum (WQS) Regression with random subset and repeated hold-out

gg <- gwqsrh(ASD ~ wqs + csex + crace + cage + YOB + OwnHome + MMC_4cat + AgeMomYrs,
             mix_name = exposures, data = finaldf, q = 10, signal = "t2",
             b = 50, rs = T,   n_vars = 12, rh = 10, validation = 0.25, 
             seed = 123123123, 
             b1_pos = T,  family = "binomial",plan_strategy = "multicore")

finaldf$wqs_q_10 <- gg$wqs

# Extracting Pearson Residuals

fit <- glm(ASD ~ wqs_q_10 + csex + crace + cage + YOB + OwnHome + MMC_4cat + AgeMomYrs, data = finaldf, family = "binomial")
finaldf$resid <-  resid(fit, type="pearson")

## SiRF algorithm with Residuals as the outcomes

n <- nrow(finaldf)

### 75 training  and 25 testing

set.seed(123456)
train.id <- sample(seq(1,n), ceiling(n*75/100))
test.id <- setdiff(1:n, train.id)

set.seed(123456)
fit <- iRF(x=finaldf[train.id, exposures], 
           y=finaldf[train.id,"resid"], 
           xtest=finaldf[test.id, exposures],
           ytest=finaldf[test.id,"resid"],
           n.iter=5, 
           n.core=1,
           select.iter = T,
           n.bootstrap=500
)

dypt <- as.data.frame(fit$interaction[fit$interaction$stability > 0.5,])
dypt$synergy <- rep(NA_real_,nrow(dypt))
for(i in 1:nrow(dypt)) dypt$synergy[i] <- sum( strsplit(dypt$int,"+")[[i]] %in% "-")
(a <- dypt[dypt$synergy == 0,])

#            int prevalence precision          cpe sta.cpe           fsd sta.fsd
# 2     DEP+_Mo+ 0.07645371 0.5584253  0.016950446   1.000  0.0055674146   0.980
# 3   DEP+_MEPB+ 0.05329762 0.5730775  0.014251962   1.000 -0.0003566124   0.464
# 5      DEP+_U+ 0.04891858 0.5474935  0.009142075   1.000  0.0003259385   0.552
# 6 DEP+_TCP246+ 0.08277110 0.4800561 -0.005390493   0.028  0.0017664920   0.684
# 7     Cd+_DEP+ 0.10020845 0.4602781 -0.015344472   0.000  0.0042578850   0.848
# 8   BPAP+_DEP+ 0.11737824 0.4493577 -0.024043531   0.000 -0.0096154233   0.008
# 
#            mip sta.mip stability synergy
# 2  0.010076128   0.912     0.936       0
# 3 -0.004036495   0.344     0.696       0
# 5  0.004338129   0.676     0.620       0
# 6 -0.017385707   0.012     0.792       0
# 7 -0.037163744   0.000     0.848       0
# 8 -0.048084104   0.000     0.928       0


### 80 training  and 20 testing

set.seed(123456)
train.id <- sample(seq(1,n), ceiling(n*80/100))
test.id <- setdiff(1:n, train.id)

set.seed(123456)
fit <- iRF(x=finaldf[train.id, exposures], 
           y=finaldf[train.id,"resid"], 
           xtest=finaldf[test.id, exposures],
           ytest=finaldf[test.id,"resid"],
           n.iter=5, 
           n.core=1,
           select.iter = T,
           n.bootstrap=500
)

dypt <- as.data.frame(fit$interaction[fit$interaction$stability > 0.5,])
dypt$synergy <- rep(NA_real_,nrow(dypt))
for(i in 1:nrow(dypt)) dypt$synergy[i] <- sum( strsplit(dypt$int,"+")[[i]] %in% "-")
(b <- dypt[dypt$synergy == 0,])

#            int prevalence precision          cpe sta.cpe           fsd sta.fsd
# 2   DEP+_MEPB+ 0.05184249 0.5652329  0.012582186   1.000 -0.0001571384   0.440
# 5 DEP+_TCP246+ 0.09818447 0.4843691 -0.004683338   0.044  0.0006587483   0.548
# 7     Cd+_DEP+ 0.11735161 0.4646551 -0.015691038   0.000  0.0058721128   0.928
# 8  Cd+_TCP246+ 0.07825674 0.4482579 -0.016512479   0.000  0.0161059772   1.000
# 
#            mip sta.mip stability synergy
# 2 -0.004100217   0.316     0.672       0
# 5 -0.018125746   0.000     0.908       0
# 7 -0.037839784   0.000     0.956       0
# 8 -0.030433645   0.000     0.512       0


### 70 training  and 30 testing

set.seed(123456)
train.id <- sample(seq(1,n), ceiling(n*70/100))
test.id <- setdiff(1:n, train.id)

set.seed(123456)
fit <- iRF(x=finaldf[train.id, exposures], 
           y=finaldf[train.id,"resid"], 
           xtest=finaldf[test.id, exposures],
           ytest=finaldf[test.id,"resid"],
           n.iter=5, 
           n.core=1,
           select.iter = T,
           n.bootstrap=500
)

dypt <- as.data.frame(fit$interaction[fit$interaction$stability > 0.5,])
dypt$synergy <- rep(NA_real_,nrow(dypt))
for(i in 1:nrow(dypt)) dypt$synergy[i] <- sum( strsplit(dypt$int,"+")[[i]] %in% "-")
(c<- dypt[dypt$synergy == 0,])

#            int prevalence precision          cpe sta.cpe          fsd sta.fsd
# 1      DEP+_U+ 0.06864853 0.5552231  0.016072350   1.000  0.001383737   0.696
# 5 DEP+_TCP246+ 0.06502159 0.4725276 -0.004330881   0.064  0.001281968   0.664
# 6     Cd+_DEP+ 0.10908710 0.4472503 -0.019703872   0.000  0.005460531   0.904
# 7   BPAP+_DEP+ 0.11549155 0.4477413 -0.020581586   0.000 -0.007251634   0.032
# 
#            mip sta.mip stability synergy
# 1  0.003823354   0.700     0.884       0
# 5 -0.014896337   0.072     0.568       0
# 6 -0.040173613   0.000     0.944       0
# 7 -0.039682690   0.000     0.928       0

intersect(intersect(a$int, b$int), c$int)
# [1] "DEP+_TCP246+" "Cd+_DEP+"


