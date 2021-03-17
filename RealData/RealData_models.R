##---------------------------------------------------------------

## Real Data Example. (Table 4,5)

## Author: Xiang Li
## Date: 03-17-2021

## Output:
## res_pq0.9.rda, res_propose_0.5.rda, res_propose_0.7.rda, res_propose_0.948.rda

##---------------------------------------------------------------




## load packages and functions
packages <- c("glmnet","MASS", "tidyverse","caret","grplasso","pROC","doParallel","foreach","plyr")
lapply(packages, require, ch=TRUE)
source("RealData_funcs.R")




## load data
load("X_Y.rda")
load("NAMCS16_Clean.RData")

## transform Y to 0, 1 (required by the grplasso)
y <- ifelse(Y==-1,0,1)

## scale the weights
X_weight <- round(dt16sub$PATWT/min(dt16sub$PATWT))

## construct lambda grid of length 125
grp_lambda_seq <- exp(seq(log(15000),0,length.out = 140))[c(1:125)]

## set seed
set.seed(2021)




## fit grplasso on data with different standardizations
cl <- parallel::makeForkCluster(2)
doParallel::registerDoParallel(cl)


stand_methods = c("none","zscore","gelman","minmax","propose")

# fit grplasso on 5 standardized data, proposed method has cutoff p=q=0.9
system.time(res_pq0.9 <- foreach(stand = stand_methods,.errorhandling = "pass",.packages = c("grplasso","pROC")) %dopar% 
              pattern_wgglasso(X, y, X_weight, grp_lambda_seq, p=0.9, q=0.9, stand))
save(res_pq0.9, file = "data/res_pq0.9.rda")

# fit grplasso on proposed standardized data with p=q=0.5
system.time(res_propose_0.5 <- foreach(stand = "propose",.errorhandling = "pass",.packages = c("grplasso","pROC")) %dopar% 
              pattern_wgglasso(X, y, X_weight, grp_lambda_seq, p=0.5, q=0.5, stand))
save(res_propose_0.5, file = "data/res_propose_0.5.rda")

# fit grplasso on proposed standardized data with p=q=0.7
system.time(res_propose_0.7 <- foreach(stand = "propose",.errorhandling = "pass",.packages = c("grplasso","pROC")) %dopar% 
              pattern_wgglasso(X, y, X_weight, grp_lambda_seq, p=0.7, q=0.7, stand))
save(res_propose_0.7, file = "data/res_propose_0.7.rda")

# fit grplasso on proposed standardized data with p=q=0.948
system.time(res_propose_0.948 <- foreach(stand = "propose",.errorhandling = "pass",.packages = c("grplasso","pROC")) %dopar% 
              pattern_wgglasso(X, y, X_weight, grp_lambda_seq, p=0.948, q=0.948, stand))
save(res_propose_0.948, file = "data/res_propose_0.948.rda")


parallel::stopCluster(cl)




## the binary prob. after dichotomizing continuous and regrouping multi-category
# proportion of binary in original data
tf_Xbi <- X[,X_bi]
summary(apply(tf_Xbi,2,function(x) max(table(x)/length(x))))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5396  0.9369  0.9818  0.9480  0.9937  0.9998 


# proposed standardization with p=q=0.5
X_dum0.5 <- lapply(X, function(x) multi_reduce(x,y=y,p=0.5,q=0.5)) %>% as.data.frame() 
summary(apply(X_dum0.5,2,function(x) max(table(x)/length(x))))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5004  0.8545  0.9704  0.8942  0.9920  0.9998 


# proposed standardization with p=q=0.7
X_dum0.7 <- lapply(X, function(x) multi_reduce(x,y=y,p=0.7,q=0.7)) %>% as.data.frame() 
summary(apply(X_dum0.7,2,function(x) max(table(x)/length(x))))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5396  0.8859  0.9708  0.9116  0.9921  0.9998 


# proposed standardization with p=q=0.9
X_dum0.9 <- lapply(X, function(x) multi_reduce(x,y=y,p=0.9,q=0.9)) %>% as.data.frame() 
summary(apply(X_dum0.9,2,function(x) max(table(x)/length(x))))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5396  0.9118  0.9714  0.9283  0.9923  0.9998 


# proposed standardization with p=q=0.948
X_dum0.948 <- lapply(X, function(x) multi_reduce(x,y=y,p=0.948,q=0.948)) %>% as.data.frame() 
summary(apply(X_dum0.948,2,function(x) max(table(x)/length(x))))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.5110  0.9280  0.9755  0.9322  0.9926  0.9998 












