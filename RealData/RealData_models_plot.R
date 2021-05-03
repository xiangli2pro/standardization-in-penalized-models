##---------------------------------------------------------------

## Plot of figure1. BIC VS. No.var & AUC  VS. No.var

## Author: Xiang Li
## Date: 03-17-2021

## Output:
## 

##---------------------------------------------------------------




## load packages and data
library(MASS)
library(tidyverse)
library(ggpubr)
library(scales)
library(ggforce)
library(patchwork)

load("data/res_pq0.9.rda")

grp_lambda_seq <- exp(seq(log(15000),0,length.out = 140))[1:125]





## organize tables
val_df <- cbind(grp_lambda_seq,t(do.call(rbind,res_pq0.9))) %>% as.data.frame()
colnames(val_df) <- c("lambda",paste0(c("logLik","AIC","BIC","AUC","No.cont","No.bi","No.multi"),
                                      "_",
                                      rep(c("none","zscore","gelman","minmax","proposed"),each=7)))

val_df <- val_df %>% 
  add_column(No.var_none = rowSums(val_df[,c("No.cont_none","No.bi_none","No.multi_none")]),
             No.var_zscore = rowSums(val_df[,c("No.cont_zscore","No.bi_zscore","No.multi_zscore")]),
             No.var_gelman = rowSums(val_df[,c("No.cont_gelman","No.bi_gelman","No.multi_gelman")]),
             No.var_minmax = rowSums(val_df[,c("No.cont_minmax","No.bi_minmax","No.multi_minmax")]),
             No.var_proposed = rowSums(val_df[,c("No.cont_proposed","No.bi_proposed","No.multi_proposed")]))

val_df_long <- val_df %>% pivot_longer(cols=-lambda,
                                       names_to=c(".value","Standardization"),
                                       names_sep = "_")

# extract tables of BIC, No.var, AUC
val_bic <- val_df %>% dplyr::select(starts_with("BIC"))
val_novar <- val_df %>% dplyr::select(starts_with("No.var"))
val_auc <- val_df %>% dplyr::select(starts_with("AUC"))

se_bic_func <- function(x){
  min_val <- min(x,na.rm = TRUE)
  se_val <- sd(x,na.rm = TRUE)
  return(which.min(abs(x-(min_val+se_val))))
}

se_auc_func <- function(x){
  max_val <- max(x,na.rm = TRUE)
  se_val <- sd(x,na.rm = TRUE)
  return(which.min(abs(x-(max_val-se_val))))
}

# bic_1se index and corresponding no.var
val_bic_midx <- apply(val_bic,2,se_bic_func)
val_bic_novar <- sapply(c(1:5),function(i) val_novar[,i][val_bic_midx[i]])
# 24 30 23 27 22

# auc_1se index and corresponding no.var
val_auc_midx <- apply(val_auc,2,se_auc_func)
val_auc_novar <- sapply(c(1:5),function(i) val_novar[,i][val_auc_midx[i]])
# 16 26 16 18 18




## plots BIC VS. No.var; AUC VS. No.var

# annotate bic_1se points of BIC
val_bic_1se_ybic <- sapply(c(1:5),function(i) val_bic[,i][val_bic_midx[i]])

p_BIC_bic1se <- ggplot(val_df_long, aes(No.var, BIC, linetype=Standardization, size=Standardization)) +
  geom_line() +
  scale_linetype_manual(breaks=c("none","zscore","gelman","minmax","proposed"),
                        values=c("longdash","solid","twodash","dashed","dotted"))+
  scale_size_manual(breaks=c("none","zscore","gelman","minmax","proposed"), 
                    values=c(0.2,0.5,0.2,0.2,0.6))+
  annotate("point", x = val_bic_novar, y = val_bic_1se_ybic, size =1.4,shape=16) +
  annotate("text", x = val_bic_novar[c(2,5)], y = val_bic_1se_ybic[c(2,5)],
           label=c("zscore","proposed"),
           vjust=c(0,0),
           hjust=c(-0.2,1.2),size=2)+
  xlim(5,80) +
  ylim(80000,100000)+
  xlab("Number of selected covariates")+
  theme_bw()+
  theme(aspect.ratio=0.8)+
  theme(legend.key.width = unit(1,"cm"))

# annotate bic_1se points of AUC
val_bic_1se_yauc <- sapply(c(1:5),function(i) val_auc[,i][val_bic_midx[i]])

p_AUC_bic1se <- ggplot(val_df_long, aes(No.var, AUC, linetype=Standardization, size=Standardization)) +
  geom_line() +
  scale_linetype_manual(breaks=c("none","zscore","gelman","minmax","proposed"),
                        values=c("longdash","solid","twodash","dashed","dotted"))+
  scale_size_manual(breaks=c("none","zscore","gelman","minmax","proposed"), 
                    values=c(0.2,0.5,0.2,0.2,0.6))+
  annotate("point", x = val_bic_novar, y = val_bic_1se_yauc, size =1.4,shape=16) +
  annotate("text", x = val_bic_novar[c(2,5)], y = val_bic_1se_yauc[c(2,5)],
           label=c("zscore","proposed"),
           vjust=c(0,0),
           hjust=c(-0.3, 1.2),size=2)+
  ylim(0.75,0.87)+
  xlim(5,80) +
  xlab("Number of selected covariates")+
  theme_bw()+
  theme(legend.key.width = unit(1,"cm"))





# annotate auc_1se points of BIC
val_auc_1se_ybic <- sapply(c(1:5),function(i) val_bic[,i][val_auc_midx[i]])

p_BIC_auc1se <- ggplot(val_df_long, aes(No.var, BIC, linetype=Standardization, size=Standardization)) +
  geom_line() +
  scale_linetype_manual(breaks=c("none","zscore","gelman","minmax","proposed"),
                        values=c("longdash","solid","twodash","dashed","dotted"))+
  scale_size_manual(breaks=c("none","zscore","gelman","minmax","proposed"), 
                    values=c(0.2,0.5,0.2,0.2,0.6))+
  annotate("point", x = val_auc_novar, y = val_auc_1se_ybic, size =1.4,shape=16) +
  annotate("text", x = val_auc_novar[c(2,5)], y = val_auc_1se_ybic[c(2,5)],
           label=c("zscore","proposed"),
           vjust=c(0,0.2),
           hjust=c(-0.2,1.2),size=2)+
  xlim(5,80) +
  xlab("Number of selected covariates")+
  ylim(80000,100000)+
  theme_bw()+
  theme(aspect.ratio=0.8)+
  theme(legend.key.width = unit(1,"cm"))


# annotate auc_1se points of AUC
val_auc_1se_yauc <- sapply(c(1:5),function(i) val_auc[,i][val_auc_midx[i]])

p_AUC_auc1se <- ggplot(val_df_long, aes(No.var, AUC, linetype=Standardization, size=Standardization)) +
  geom_line() +
  scale_linetype_manual(breaks=c("none","zscore","gelman","minmax","proposed"),
                        values=c("longdash","solid","twodash","dashed","dotted"))+
  scale_size_manual(breaks=c("none","zscore","gelman","minmax","proposed"), 
                    values=c(0.2,0.5,0.2,0.2,0.6))+
  annotate("point", x = val_auc_novar, y = val_auc_1se_yauc, size =1.4,shape=16) +
  annotate("text", x = val_auc_novar[c(2,5)], y = val_auc_1se_yauc[c(2,5)],
           label=c("zscore","proposed"),
           vjust=c(0,0),
           hjust=c(-0.3, 1.1),size=2)+
  ylim(0.75,0.87)+
  xlim(5,80) +
  xlab("Number of selected covariates")+
  theme_bw()+
  theme(legend.key.width = unit(1,"cm"))



