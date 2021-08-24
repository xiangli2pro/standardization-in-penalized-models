##---------------------------------------------------------------

## Conduct Simulation 1. (Table 1)

## Author: Xiang Li
## Date: 05-13-2021

## Output:
## res_coef_update.rda (simulation raw table)
## res_coef_t2.rda (summary mean(sd) coefficients)
## res_ratio_t2.rda (summary ratio(sd) coefficients)

##---------------------------------------------------------------




## load packages and functions
packages <- c("glmnet","MASS", "tidyverse","doParallel","foreach","boot","faraway","caret")
lapply(packages, require, ch=TRUE)
source("Simulation_funcs.R")




## parameters of the normal distribution of data X
n <- 10000 # number of observations
pn <- 3 # number of covariates
mu <- rep(0,pn)
Sigma <- diag(rep(1,pn))
set.seed(2020)

## true coefficients
beta_0 <- 4
beta_p <- c(2,2,-2)





## function: return standardized coefficients under different standardization methods
coeff_func_p <- function(beta_0, beta_p, n, mu,Sigma, p=0.5){
  
  # generate latent data X
  X <- mvrnorm(n, mu=mu, Sigma = Sigma) %>% as.data.frame()
  colnames(X) <- paste0("X",c(1,2,3))
  
  # true model
  linear_part <- beta_0 + (as.matrix(X) %*% beta_p)
  Y <- rbinom(nrow(X),size = 1, prob = ilogit(linear_part))
  
  # generate observed data (dichotomize X2 to Z2)
  Z <- data.frame(Z1=X$X1,
                  Z2=ifelse(X$X2<=qnorm(p,mu[2],sd = sqrt(Sigma[2,2])),0,1),
                  Z3=factor(ifelse(X$X3<=qnorm(1/3,mu[3],sd = sqrt(Sigma[3,3])),"f0",
                                   ifelse(X$X3<=qnorm(2/3,mu[3],sd = sqrt(Sigma[3,3])),"f1","f2")),levels = paste0("f",c(0,1,2)))
  )
  
  # encode Z
  dummies <- dummyVars(~., data = Z, fullRank = T)
  Z_dum <- predict(dummies,Z)
  ZY_data <- cbind(Y,Z_dum) 
  ZY_data_prior <- cbind(Y,Z)
  
  # glm model on data with no standardization
  lg_mod1 <- ZY_data %>% as.data.frame() %>% glm(Y~., data = ., family = "binomial")
  
  # glm model on data with Z-score standardization
  lg_mod2 <- ZY_data %>% as.data.frame() %>% mutate_at(c(2:5),function(x) x/sd(x)) %>% glm(Y~., data = ., family = "binomial")
  
  # glm model on data with Gelman standardization
  lg_mod3 <- ZY_data %>% as.data.frame() %>% mutate_at(c(2),function(x) x/2/sd(x)) %>% glm(Y~., data = ., family = "binomial")
  
  # glm model on data with Min-max standardization
  lg_mod4 <- ZY_data %>% as.data.frame() %>% mutate_at(c(2),function(x) (x-min(x))/(max(x)-min(x))) %>% glm(Y~., data = ., family = "binomial")
  
  # glm model on data with Proposed standardization, no 2sd scale
  lg_mod5 <- ZY_data_prior %>% 
    mutate(Z1=ifelse(Z1<=quantile(Z1,probs = mean(Z2==0)),0,1)) %>%  
    mutate_at(c(4),function(x) multi_fix(x)) %>% 
    glm(Y~., data = ., family = "binomial")
  
  # glm model on data with Proposed standardization, with 2sd scale
  lg_mod6 <- ZY_data_prior %>% 
    mutate(Z1=ifelse(Z1<=quantile(Z1,probs = mean(Z2==0)),0,1)) %>%  
    mutate_at(c(4),function(x) multi_fix(x)) %>% 
    mutate_at(c(2:4), function(x) return(x/2/sd(x))) %>% 
    glm(Y~., data = ., family = "binomial")
  
  # get coefficients from all models
  coef_stat <- lapply(paste0("lg_mod",c(1:6)), function(x) summary(get(x))$coefficients[,1]) %>% unlist()
  
  return(c(coef_stat))
}





## Run 2000 simulations
## open parallel computing clusters 
cl <- makeForkCluster(2)
registerDoParallel(cl)

sim_num <- 2000
# data when X_bi is dichotomized by p=0.5
system.time(res_p0.5 <- foreach(m=c(1:sim_num),.combine = "rbind",.packages = c("glmnet","dplyr", "MASS")) %dopar% 
              coeff_func_p(beta_0, beta_p,n, mu,Sigma, p=0.5))
# data when X_bi is dichotomized by p=0.7
system.time(res_p0.7 <- foreach(m=c(1:sim_num),.combine = "rbind",.packages = c("glmnet","dplyr", "MASS")) %dopar% 
              coeff_func_p(beta_0, beta_p,n, mu,Sigma, p=0.7))
# data when X_bi is dichotomized by p=0.9
system.time(res_p0.9 <- foreach(m=c(1:sim_num),.combine = "rbind",.packages = c("glmnet","dplyr", "MASS")) %dopar% 
              coeff_func_p(beta_0, beta_p,n, mu,Sigma, p=0.9))

# save data
save(res_p0.5, res_p0.7,res_p0.9, file = "simulation1_data/res_coef_update.rda")




## summary (mean, sd, ratio) table 
## outlier percentile 0.95
q <- 0.95

## mean(sd) table
mean_table <- cbind(sum_func(res_p0.5,q)$p_mean,
                    sum_func(res_p0.7,q)$p_mean,
                    sum_func(res_p0.9,q)$p_mean) %>% round(2)

sd_table <- cbind(sum_func(res_p0.5,q)$p_sd,
                  sum_func(res_p0.7,q)$p_sd,
                  sum_func(res_p0.9,q)$p_sd) %>% round(2)

# paste the mean (sd)
res_coef_t1 <- matrix(paste0(mean_table," (",sd_table,")"),nrow = nrow(mean_table),byrow = FALSE)

# create row names
row_name <- c(rep(c("Intercept","Z_cont","Z_bi","Z_multi_f1","Z_multi_f2"),4),
              rep(c("Intercept","Z_cont","Z_bi","Z_multi"),2))

# add coeff. names as a new column
res_coef_t2 <- cbind(Coeff=row_name,res_coef_t1)
colnames(res_coef_t2) <- c("Coeff","p=0.5","p=0.7","p=0.9")

# save data
save(res_coef_t2,file = "simulation1_data/res_coef_t2.rda")

## ratio(sd) table
ratio_mean_table <- cbind(sum_func(res_p0.5,q)$p_ratio_mean,
                          sum_func(res_p0.7,q)$p_ratio_mean,
                          sum_func(res_p0.9,q)$p_ratio_mean) %>% round(2)

ratio_sd_table <- cbind(sum_func(res_p0.5,q)$p_ratio_sd,
                        sum_func(res_p0.7,q)$p_ratio_sd,
                        sum_func(res_p0.9,q)$p_ratio_sd) %>% round(2)

ratio_CM_mean_table <- cbind(sum_func(res_p0.5,q)$pCM_ratio_mean,
                             sum_func(res_p0.7,q)$pCM_ratio_mean,
                             sum_func(res_p0.9,q)$pCM_ratio_mean) %>% round(2)

ratio_CM_sd_table <- cbind(sum_func(res_p0.5,q)$pCM_ratio_sd,
                           sum_func(res_p0.7,q)$pCM_ratio_sd,
                           sum_func(res_p0.9,q)$pCM_ratio_sd) %>% round(2)

ratio_BM_mean_table <- cbind(sum_func(res_p0.5,q)$pBM_ratio_mean,
                             sum_func(res_p0.7,q)$pBM_ratio_mean,
                             sum_func(res_p0.9,q)$pBM_ratio_mean) %>% round(2)

ratio_BM_sd_table <- cbind(sum_func(res_p0.5,q)$pBM_ratio_sd,
                           sum_func(res_p0.7,q)$pBM_ratio_sd,
                           sum_func(res_p0.9,q)$pBM_ratio_sd) %>% round(2)

res_ratio_t1 <- matrix(paste0(ratio_mean_table," (",ratio_sd_table,")"),nrow = nrow(ratio_mean_table),byrow = FALSE)
res_ratio_CM <- matrix(paste0(ratio_CM_mean_table," (",ratio_CM_sd_table,")"),nrow = nrow(ratio_CM_mean_table),byrow = FALSE)
res_ratio_BM <- matrix(paste0(ratio_BM_mean_table," (",ratio_BM_sd_table,")"),nrow = nrow(ratio_BM_mean_table),byrow = FALSE)

res_ratio_t2 <- rbind(res_ratio_t1,res_ratio_CM, res_ratio_BM)
res_ratio_t2 <- cbind(c(rep("Z_cont/Z_bi", 6),rep("Z_cont/Z_multi",2),rep("Z_bi/Z_multi",2)),
                      res_ratio_t2)
colnames(res_ratio_t2) <- c("Coeff_ratio","p=0.5","p=0.7","p=0.9")

# save data
save(res_ratio_t2,file = "simulation1_data/res_ratio_t2.rda")

