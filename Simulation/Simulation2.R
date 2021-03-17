##---------------------------------------------------------------

## Conduct Simulation 2. (Table 2,3)

## Author: Xiang Li
## Date: 03-16-2021

## Output:
##


##---------------------------------------------------------------




## load packages and functions
packages <- c("glmnet","MASS", "tidyverse","doParallel","foreach","caret","faraway","gglasso")
lapply(packages, require, ch=TRUE)
source("Simulation_funcs.R")




## true coefficients
beta_0 <- 0
beta_p <- c(rep(0.05,3),rep(-0.05,3),rep(0,6))




##----------------------
## functions
##----------------------

## functions: return standardized coefficients under different standardization methods
penal_logist <- function(beta_0, beta_p,K=3, n,p,mu,Sigma, dist="normal"){
  
  if (dist=="normal"){
    
    # generate latent data X
    X <- mvrnorm(n = n,mu = mu,Sigma = Sigma) %>% as.data.frame()
    colnames(X) <- paste0("X",c(1:p))
    
    # transform to observed data
    Z <- Z_gen_normal(X)
    
  } else if (dist=="exp"){
    
    rate <- rep(1,p)
    X <- sapply(rate, function(x) rexp(n,x)) %>% as.data.frame()
    colnames(X) <- paste0("X",c(1:p))
    
    Z <- Z_gen_exp(X)
    
  } else if (dist=="gamma"){
    
    g_alpha <- 10
    g_beta <- 0.32
    X <- sapply(c(1:p), function(x) rgamma(n,shape=g_alpha,scale=g_beta)) %>% as.data.frame()
    colnames(X) <- paste0("X",c(1:p))
    
    Z <- Z_gen_gamma(X)
    
  }
  
  # true model
  linear_part <- beta_0 + (as.matrix(X) %*% beta_p)
  Y <- rbinom(nrow(X),size = 1, prob = ilogit(linear_part))
  Y_np1 <- ifelse(Y==0,-1,1)

  # group index
  group_z <- c(c(1:4),rep(5,2),rep(6,3),c(7:10),rep(11,2),rep(12,3))
  
  # encode Z
  dummies1 <- dummyVars(~., data = Z, fullRank = T)
  Z_dum1 <- predict(dummies1,Z)
  
  # 1. gglasso model on data with no standardization
  mod_1 <- cv.gglasso(x=Z_dum1,y=Y_np1,group=group_z,loss="logit",pred.loss = "misclass",nfolds = K)
  
  # 2. gglasso model on data with Gelman standardization
  Z_dum2 <- Z_dum1 %>% as.data.frame() %>% 
    mutate_at(c('Z1','Z7'),function(x) return(x/2/sd(x))) %>% 
    as.matrix()
  
  mod_2 <- cv.gglasso(x=Z_dum2,y=Y_np1,group=group_z,loss="logit",pred.loss = "misclass",nfolds = K)
  
  # 3. gglasso model on data with Min-max standardization
  Z_dum3 <- Z_dum1 %>% as.data.frame() %>% 
    mutate_at(c('Z1','Z7'),function(x) return((x-min(x))/(max(x)-min(x)))) %>% 
    as.matrix()
  
  mod_3 <- cv.gglasso(x=Z_dum3,y=Y_np1,group=group_z,loss="logit",pred.loss = "misclass",nfolds = K)
  
  # 4. gglasso model on data with Proposed standardization, with no 2sd
  Z_bi4 <- Z %>% 
    mutate(Z1=ifelse(Z1<=quantile(Z1,probs = 0.5),0,1)) %>% 
    mutate(Z7=ifelse(Z7<=quantile(Z7,probs = 0.5),0,1)) %>% 
    mutate_at(c('Z5','Z6'), function(x) multi_fix(x)) %>%
    mutate_at(c('Z11','Z12'), function(x) multi_random(x)) %>%as.matrix()
  
  mod_4 <- cv.gglasso(x=Z_bi4,y=Y_np1,group=c(1:ncol(X)),loss="logit",pred.loss = "misclass",nfolds = K)
  
  # 5. gglasso model on data with Proposed standardization, with 2sd
  Z_bi5 <- Z_bi4 %>% as.data.frame() %>%  
    mutate_all(., function(x) return(x/2/sd(x))) %>% 
    as.matrix()
  
  mod_5 <- cv.gglasso(x=Z_bi5,y=Y_np1,group=c(1:ncol(X)),loss="logit",pred.loss = "misclass",nfolds = K)
  
  # 6. gglasso model on data with Z-score standardization
  Z_dum6 <- Z_dum1 %>% as.data.frame() %>% 
    mutate_all(.,function(x) return(x/sd(x))) %>% 
    as.matrix()
  
  mod_6 <- cv.gglasso(x=Z_dum6,y=Y_np1,group=group_z,loss="logit",pred.loss = "misclass",nfolds = K)
  
  # summary
  mod_idx <- c(1:6)
  lambda_1se <- sapply(paste0("mod_",mod_idx), function(x) get(x)$lambda.1se) ## 6
  sens_1se <- sapply(paste0("mod_",mod_idx), function(x) get_sens_1se(get(x),group_z)) 
  spec_1se <- sapply(paste0("mod_",mod_idx), function(x) get_spec_1se(get(x),group_z)) 
  cvm_1se <- sapply(paste0("mod_",mod_idx), function(x) get_cvm_1se(get(x))) 
  
  lambda_min <- sapply(paste0("mod_",mod_idx), function(x) get(x)$lambda.min)
  sens_min <- sapply(paste0("mod_",mod_idx), function(x) get_sens_min(get(x),group_z))
  spec_min <- sapply(paste0("mod_",mod_idx), function(x) get_spec_min(get(x),group_z))
  cvm_min <- sapply(paste0("mod_",mod_idx), function(x) get_cvm_min(get(x)))
  
  select_1se <- c(sapply(paste0("mod_",mod_idx), function(x) get_select_1se(get(x),group_z))) ## 72
  select_min <- c(sapply(paste0("mod_",mod_idx), function(x) get_select_min(get(x),group_z))) 
  
  return(c(lambda_1se,sens_1se,spec_1se,cvm_1se,
           lambda_min,sens_min,spec_min,cvm_min,
           select_1se,select_min))
}




##----------------------
## 2000 simulations
##----------------------

## parameters of normal distribution: independent
n <- 10000
p <- 12
mu <- rep(0,p)
Sigma <- diag(rep(1,p))
set.seed(2021)

system.time(
  res_norm_indep <- replicate(n=2000, penal_logist(beta_0, beta_p,K=3,n,p,mu,Sigma, dist="normal"), simplify = FALSE)
)
save(res_norm_indep,file = "simulation2_data/res_norm_indep.rda")




## parameters of normal distribution: partially correlated
Sigma <- diag(rep(1,p))
rho <- 0.3
for (i in c(1:6)){
  for (j in c(1:6)){
    if (i!=j) Sigma[i,j] = rho
  }
}

system.time(
  res_norm_part <- replicate(n=2000, penal_logist(beta_0, beta_p,K=3,n,p,mu,Sigma, dist="normal"), simplify = FALSE)
)
save(res_norm_part,file = "simulation2_data/res_norm_part.rda")




## parameters of normal distribution: all correlated
Sigma <- diag(rep(1,p))
rho <- 0.3
for (i in c(1:p)){
  for (j in c(1:p)){
    if (i!=j) Sigma[i,j] = rho
  }
}

system.time(
  res_norm_all <- replicate(n=2000, penal_logist(beta_0, beta_p,K=3,n,p,mu,Sigma, dist="normal"), simplify = FALSE)
)
save(res_norm_all,file = "simulation2_data/res_norm_all.rda")




## latent exponential distribution
system.time(
  res_exp <- replicate(n=2000, penal_logist(beta_0, beta_p,K=3,n,p,mu,Sigma, dist="exp"), simplify = FALSE)
)
save(res_exp,file = "simulation2_data/res_exp.rda")




## latent gamma distribution
system.time(
  res_gamma <- replicate(n=2000, penal_logist(beta_0, beta_p,K=3,n,p,mu,Sigma, dist="gamma"), simplify = FALSE)
)
save(res_gamma,file = "simulation2_data/res_gamma.rda")




##----------------------
## summary tables
##----------------------
tb_norm_indep <- summary_func(res_norm_indep)
tb_norm_part <- summary_func(res_norm_part)
tb_norm_all <- summary_func(res_norm_all)
tb_exp <- summary_func(res_exp)
tb_gamma <- summary_func(res_gamma)

save(tb_norm_indep, tb_norm_part, tb_norm_all, tb_exp, tb_gamma, file="simulation2_data/summary_tables.rda")




