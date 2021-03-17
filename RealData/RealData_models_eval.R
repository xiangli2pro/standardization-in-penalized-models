##---------------------------------------------------------------

## Evaluate performance (AIC, BIC, AUC, etc) of models generated in RealData_models.R

## Author: Xiang Li
## Date: 03-17-2021

## Output:
## 

##---------------------------------------------------------------




## load packages and functions
packages <- c("glmnet","MASS", "tidyverse","caret","grplasso","pROC","doParallel","foreach","plyr")
lapply(packages, require, ch=TRUE)
source("RealData_funcs.R")

grp_lambda_seq <- exp(seq(log(15000),0,length.out = 140))[c(1:125)]



##----------------------
## get AIC, BIC, logLik, AUC of models in res_pq0.9.rda
##----------------------
load("data/res_pq0.9.rda")

# turn list to data.frame
val_df <- cbind(grp_lambda_seq,t(do.call(rbind,res_pq0.9))) %>% as.data.frame()
colnames(val_df) <- c("lambda",paste0(c("logLik","AIC","BIC","AUC",
                                        "No.cont","No.bi","No.multi"),
                                      "_",
                                      rep(c("none","zscore","gelman","minmax","propose"),each=7)))

# add column number of total selected covariates for each standardization
val_df <- val_df %>% 
  add_column(No.var_none = rowSums(val_df[,c("No.cont_none","No.bi_none","No.multi_none")]),
             No.var_zscore = rowSums(val_df[,c("No.cont_zscore","No.bi_zscore","No.multi_zscore")]),
             No.var_gelman = rowSums(val_df[,c("No.cont_gelman","No.bi_gelman","No.multi_gelman")]),
             No.var_minmax = rowSums(val_df[,c("No.cont_minmax","No.bi_minmax","No.multi_minmax")]),
             No.var_propose = rowSums(val_df[,c("No.cont_propose","No.bi_propose","No.multi_propose")]))

# extrac data of BIC, No.var, AUC and logLik
val_bic <- val_df %>% dplyr::select(starts_with("BIC"))
val_novar <- val_df %>% dplyr::select(starts_with("No.var"))
val_auc <- val_df %>% dplyr::select(starts_with("AUC"))
val_logLik <- val_df %>% dplyr::select(starts_with("logLik"))




## informations corresponding to bic_1se
# for each stand. find the corresponding bic_1se index
val_bic_midx <- apply(val_bic,2,se_bic_func)

# for each stand. find the --novar-- corresponding to their bic_1se 
val_bic_novar <- sapply(c(1:5),function(i) val_novar[,i][val_bic_midx[i]])
# 24 30 23 27 22

# # for each stand. find the --lambda-- corresponding to their bic_1se 
# grp_lambda_seq[val_bic_midx]
# # 383.5029 1529.8369  383.5029  357.8697 + 1160.0317 1010.1415 1082.4953
# 
# # for each stand. find the --bic-- corresponding to their bic_1se 
# round(sapply(c(1:5),function(i) val_bic[,i][val_bic_midx[i]]),3)
# # 90189.60 85915.15 90144.88 86735.93 + 85290.75 85120.22 85065.39
# 
# # for each stand. find the --auc-- corresponding to their bic_1se 
# sapply(c(1:5),function(i) val_auc[,i][val_bic_midx[i]])
# # 0.8065029 0.8360155 0.8071335 0.8316302 + 0.8405681 0.8414804 0.8415315
# 
# # for each stand. find the --logLik-- corresponding to their bic_1se 
# sapply(c(1:5),function(i) val_logLik[,i][val_bic_midx[i]])
# # -44896.71 -42699.60 -44860.53 -43059.32 + -42345.94 -42205.39 -42251.68




## informations corresponding to auc_1se
val_auc_midx <- apply(val_auc,2,se_auc_func)
val_auc_novar <- sapply(c(1:5),function(i) val_novar[,i][val_auc_midx[i]])
# 16 26 16 18 18





##----------------------
## use the lambda corresponding to bic_1se and auc_1se to refit grplasso
## extract the covariate names
##----------------------
load("X_Y.rda")
load("NAMCS16_Clean.RData")
y <- ifelse(Y==-1,0,1)
X_weight <- round(dt16sub$PATWT/min(dt16sub$PATWT))

##refit grplasso 
cl <- parallel::makeForkCluster(2)
doParallel::registerDoParallel(cl)

stand_methods = c("none","zscore","gelman","minmax","propose")

system.time(features_bic <- foreach(m=1:5,.errorhandling = "pass",.packages = c("grplasso","pROC")) %dopar% 
              get_features(X, y, X_weight,
                           lambda=grp_lambda_seq[val_bic_midx[m]],
                           p=0.9, q=0.9, stand=stand_methods[m]))

system.time(features_auc <- foreach(m=1:5,.errorhandling = "pass",.packages = c("grplasso","pROC")) %dopar% 
              get_features(X, y, X_weight,
                           lambda=grp_lambda_seq[val_auc_midx[m]],
                           p=0.9, q=0.9, stand=stand_methods[m]))

parallel::stopCluster(cl)




##----------------------
## inspect the selected covariates by proposed method when p=q=0.5, p=q=0.7, p=q=0.948
## use the lambda corresponding to bic_1se and auc_1se to refit grplasso
## extract the covariate names
##----------------------
load("data/res_propose_0.5.rda")
load("data/res_propose_0.7.rda")
load("data/res_propose_0.948.rda")

prop_features <- function(prop_df, p=0.5){
  
  val_prop_df <- cbind(grp_lambda_seq, t(prop_df[[1]])) %>% as.data.frame()
  colnames(val_prop_df) <- c("lambda",paste0(c("logLik","AIC","BIC","AUC","No.cont","No.bi","No.multi"),
                                             "_",
                                             "prop"))
  val_prop_df <- val_prop_df %>% 
    add_column(No.var_prop = rowSums(val_prop_df[, c("No.cont_prop","No.bi_prop","No.multi_prop")]))
  
  val_bic_midx <- se_bic_func(val_prop_df$BIC_prop)
  val_bic_novar <- val_prop_df$No.var_prop[val_bic_midx]
  
  val_auc_midx <- se_auc_func(val_prop_df$AUC_prop)
  val_auc_novar <- val_prop_df$No.var_prop[val_auc_midx]
  
  features_bic <- get_features(X, y, X_weight,
                               lambda=grp_lambda_seq[val_bic_midx],p=p, q=p, stand="propose")
  features_auc <- get_features(X, y, X_weight,
                               lambda=grp_lambda_seq[val_auc_midx],p=p, q=p, stand="propose")
  
  return(list(features_bic=features_bic, features_auc=features_auc,
              val_bic_midx=val_bic_midx, val_auc_midx=val_auc_midx))
}

features_p0.5 <- prop_features(res_propose_0.5, p=0.5)
features_p0.7 <- prop_features(res_propose_0.7, p=0.7)
features_p0.948 <- prop_features(res_propose_0.948, p=0.948)

library(rlist)
features_bic <- list.append(features_bic, 
                            features_p0.5$features_bic,
                            features_p0.7$features_bic,
                            features_p0.948$features_bic)

features_auc <- list.append(features_auc, 
                            features_p0.5$features_auc,
                            features_p0.7$features_auc,
                            features_p0.948$features_auc)

val_bic_midx <- c(val_bic_midx, 
                  features_p0.5$val_bic_midx,
                  features_p0.7$val_bic_midx,
                  features_p0.948$val_bic_midx)

val_auc_midx <- c(val_auc_midx, 
                  features_p0.5$val_auc_midx,
                  features_p0.7$val_auc_midx,
                  features_p0.948$val_auc_midx)

save(features_bic, features_auc, val_bic_midx, val_auc_midx, file = "data/res_features.rda")

tf <- do.call(cbind.fill, features_bic)
write.csv(tf, file='features_bic.xlsx')




## No.cont, No.bi, No.multi
## 0.5 
c(sum(features_bic[[6]] %in% X_cont),sum(features_bic[[6]] %in% X_bi),sum(features_bic[[6]] %in% X_multi))
# 1 17  9

## 0.7
c(sum(features_bic[[7]] %in% X_cont),sum(features_bic[[7]] %in% X_bi),sum(features_bic[[7]] %in% X_multi))
# 1 10  7

## 0.9
c(sum(features_bic[[5]] %in% X_cont),sum(features_bic[[5]] %in% X_bi),sum(features_bic[[5]] %in% X_multi))
# 1 12  9
