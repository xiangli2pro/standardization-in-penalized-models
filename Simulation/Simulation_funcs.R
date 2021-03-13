##---------------------------------------------------------------

## Functions called in Simulation1.R and Simulation2.R

## Author: Xiang Li
## Date: 03-12-2021

## Functions:
## sum_func(), multi_fix(), multi_random()
## nzero_coef_1se(), get_select_1se(), get_sens_1se(), get_sepc_1se(), get_cvm_1se()
## nzero_coef_min(), get_select_min(), get_sens_min(), get_sepc_min(), get_cvm_min()
## mod_plot(), plot_data(), Z_gen(), Z_exp_gamma(), Z_gen_gamma()

##---------------------------------------------------------------




## summary table after removing outliers beyond 95th percent values
sum_func <- function(p_data, q=0.95){
  
  rm_outlier <- function(x,q) {
      x %>% as.data.frame() %>% 
      filter(Z1<=quantile(Z1,probs = q)) %>% 
      filter(Z2<=quantile(Z2,probs = q))
  }
  
  mean_func <- function(x,q){
    x <- rm_outlier(x,q)
    x_mean <- colMeans(x)
    return(c(x_mean))
  }
  
  sd_func <- function(x,q){
    x <- rm_outlier(x,q)
    x_sd <- apply(x, 2, sd)
    return(c(x_sd))
  }
  
  # ratio of Z_cont / Z_bi
  ratio_mean <- function(x,q){
    x <- rm_outlier(x,q)
    x = as.data.frame(x) 
    return(mean(x$Z1/x$Z2))
  }
  
  ratio_sd <- function(x,q){
    x <- rm_outlier(x,q)
    x = as.data.frame(x) 
    return(sd(x$Z1/x$Z2))
  }
  
  # ratio of Z_cont / Z_multi
  ratio_CM_mean <- function(x,q){
    #x_data <- rm_outlier(x,q)
    x = as.data.frame(x) 
    return(mean(x$Z1/x$Z3))
  }
  
  ratio_CM_sd <- function(x,q){
    #x_data <- rm_outlier(x,q)
    x = as.data.frame(x) 
    return(sd(x$Z1/x$Z3))
  }
  
  # ratio of Z_bi / Z_multi
  ratio_BM_mean <- function(x,q){
    #x_data <- rm_outlier(x,q)
    x = as.data.frame(x) 
    return(mean(x$Z2/x$Z3))
  }
  
  ratio_BM_sd <- function(x,q){
    #x_data <- rm_outlier(x,q)
    x = as.data.frame(x) 
    return(sd(x$Z2/x$Z3))
  }
  
  # mean coefficients of data standardized by different approaches
  p_mean <- unlist(sapply(list(p_data[,c(1:5)], p_data[,c(6:10)], 
                               p_data[,c(11:15)], p_data[,c(16:20)],
                               p_data[,c(21:24)],p_data[,c(25:28)]),
                          function(x) mean_func(x,q)))
  
  # sd coefficients of data standardized by different approaches
  p_sd <- unlist(sapply(list(p_data[,c(1:5)], p_data[,c(6:10)], 
                             p_data[,c(11:15)], p_data[,c(16:20)],
                             p_data[,c(21:24)],p_data[,c(25:28)]),
                        function(x) sd_func(x,q)))
  
  # mean coefficients ratio of Z_cont/Z_bi of data standardized by different approaches
  p_ratio_mean <- unlist(sapply(list(p_data[,c(1:5)], p_data[,c(6:10)], 
                                     p_data[,c(11:15)], p_data[,c(16:20)],
                                     p_data[,c(21:24)],p_data[,c(25:28)]),
                                function(x) ratio_mean(x,q)))
  
  p_ratio_sd <- unlist(sapply(list(p_data[,c(1:5)], p_data[,c(6:10)], 
                                   p_data[,c(11:15)], p_data[,c(16:20)],
                                   p_data[,c(21:24)],p_data[,c(25:28)]),
                              function(x) ratio_sd(x,q)))
  
  # mean coefficients ratio of Z_cont/Z_multi of data standardized by proposed method
  pCM_ratio_mean <- unlist(sapply(list(p_data[,c(21:24)],p_data[,c(25:28)]),
                                  function(x) ratio_CM_mean(x,q)))
  
  pCM_ratio_sd <- unlist(sapply(list(p_data[,c(21:24)],p_data[,c(25:28)]),
                                function(x) ratio_CM_sd(x,q)))
  
  # mean coefficients ratio of Z_bi/Z_multi of data standardized by proposed method
  pBM_ratio_mean <- unlist(sapply(list(p_data[,c(21:24)],p_data[,c(25:28)]),
                                  function(x) ratio_BM_mean(x,q)))
  
  pBM_ratio_sd <- unlist(sapply(list(p_data[,c(21:24)],p_data[,c(25:28)]),
                                function(x) ratio_BM_sd(x,q)))
  
  # return results
  return(list(p_mean=p_mean,p_sd=p_sd,
              p_ratio_mean=p_ratio_mean,p_ratio_sd=p_ratio_sd,
              pCM_ratio_mean=pCM_ratio_mean, pCM_ratio_sd=pCM_ratio_sd,
              pBM_ratio_mean=pBM_ratio_mean, pBM_ratio_sd=pBM_ratio_sd))
}





## Proposed standardization: reduce true multi-factor to binary 
## by regrouping the fixed levels in ascending order
multi_fix <- function(x){
  
  # for multi-factor with 3 levels
  if (length(levels(x))==3){
    x_levels <- levels(x)
    g0_index <- x %in% (x_levels)[c(1)]
    
    new_x <- vector(length = length(x))
    new_x[g0_index] <- 1 
    new_x[!g0_index] <- 0 
    return(new_x)
    
  # for multi-factor with 4 levels
  }else if (length(levels(x))==4){
    x_levels <- levels(x)
    g0_index <- x %in% (x_levels)[c(1,2)]
    
    new_x <- vector(length = length(x))
    new_x[g0_index] <- 1
    new_x[!g0_index] <- 0 
    return(new_x)
  } else {
    
    return(x)
  }
}




## Proposed standardization: reduce noise multi-factor to binary 
## by regrouping the fixed levels randomly
multi_random <- function(x){

  # for multi-factor with 3 levels
  if (length(levels(x))==3){
    r_group <- sample(c(1:3),size = 3, replace = FALSE)
    x_levels <- levels(x)
    
    g0_index <- x %in% (x_levels)[r_group[1]]
    new_x <- vector(length = length(x))
    new_x[g0_index] <- 1
    new_x[!g0_index] <- 0 
    
    return(new_x)
    
  # for multi-factor with 4 levels
  }else if (length(levels(x))==4){
    
    r_group <- sample(c(1:4),size = 4, replace = FALSE)
    x_levels <- levels(x)
    
    g0_index <- x %in% (x_levels)[r_group[c(1,2)]]
    new_x <- vector(length = length(x))
    new_x[g0_index] <- 1
    new_x[!g0_index] <- 0
    
    return(new_x)
    
  } else {
    return(x)
  }
}





## cv = cv_1se: return non-zero coefficient index
nzero_coef_1se <- function(mod, group_z){
  
  nzero_coef <- which(coef(mod$gglasso.fit, s=mod$lambda.1se)[-1]!=0)
  
  if (length(mod$gglasso.fit$group)==18){
    nzero_ind <- unique(group_z[nzero_coef])
  } else {
    nzero_ind <- nzero_coef
  }
  
  return(nzero_ind)
}




## cv = cv_1se: return selection indicator 1:selected, 0:not selected
get_select_1se <- function(mod,group_z){
  
  nzero_ind <- nzero_coef_1se(mod,group_z)
  
  return(c(1:12) %in% nzero_ind)
}




## cv = cv_1se: return sensitivity of true covariates
get_sens_1se <- function(mod,group_z){
  
  nzero_ind <- nzero_coef_1se(mod,group_z)
  
  n_sens <- sum(c(1:6) %in% nzero_ind)/6
  return(n_sens)
}




## cv = cv_1se:  return sepcificity of noise covariates
get_spec_1se <- function(mod,group_z){
  
  nzero_ind <- nzero_coef_1se(mod,group_z)
  
  # rate of not being selected for noise
  n_spec <- sum(!c(7:12) %in% nzero_ind)/6
  return(n_spec)
}




## cv = cv_1se: return cv
get_cvm_1se <- function(mod){
  mod$cvm[which(mod$lambda==mod$lambda.1se)]
}




## cv = cv_min: return non-zero coefficient index
nzero_coef_min <- function(mod,group_z){
  
  nzero_coef <- which(coef(mod$gglasso.fit, s=mod$lambda.min)[-1]!=0)
  
  if (length(mod$gglasso.fit$group)==18){
    nzero_ind <- unique(group_z[nzero_coef])
  } else {
    nzero_ind <- nzero_coef
  }
  
  return(nzero_ind)
}




## cv = cv_min: return selection index
get_select_min <- function(mod,group_z){
  
  nzero_ind <- nzero_coef_min(mod,group_z)
  
  return(c(1:12) %in% nzero_ind)
}




## cv = cv_min: return sensitivity
get_sens_min <- function(mod,group_z){
  
  nzero_ind <- nzero_coef_min(mod,group_z)
  
  n_sens <- sum(c(1:6) %in% nzero_ind)/6
  return(n_sens)
}




## cv = cv_min: return specificity
get_spec_min <- function(mod,group_z){
  
  nzero_ind <- nzero_coef_min(mod,group_z)
  
  n_spec <- sum(!c(7:12) %in% nzero_ind)/6
  return(n_spec)
}




## cv = cv_min: return cv_min
get_cvm_min <- function(mod){
  mod$cvm[which(mod$lambda==mod$lambda.min)]
}




## plot solution path of coefficients from gglasso model
mod_plot <- function(mod, mod_title, mod_select="None"){
  
  # get coefficients for all lambda candidates
  beta <- mod$gglasso.fit$beta
  coef_num <- ifelse(nrow(beta)==18,9,6)
  
  # rearrange to long format
  mod_coef <- as.data.frame(beta)
  mod_coef$coef <- coef_names <-  row.names(mod_coef)
  mod_coef$true <- c(rep("true",coef_num),rep("noise",coef_num))
  mod_coef <- reshape::melt(mod_coef, id=c("coef","true"))
  mod_coef$variable <- as.numeric(gsub("s","",mod_coef$variable))
  mod_coef$lambda <- mod$lambda[mod_coef$variable+1]
  mod_coef$norm <- apply(abs(beta),2,sum)[mod_coef$variable+1]
  mod_df_1se <- mod$gglasso.fit$df[which(mod$lambda==mod$lambda.1se)]
  mod_df_min <- mod$gglasso.fit$df[which(mod$lambda==mod$lambda.min)]
  
  # select true coefficients or noise coefficients
  if(mod_select=="true"){
    mod_coef <- mod_coef %>% filter(coef %in% coef_names[1:coef_num])
  } else if (mod_select=="noise"){
    mod_coef <- mod_coef %>% filter(! coef %in% coef_names[1:coef_num])
  }
  
  ggplot(mod_coef, aes(lambda, value,color=coef,linetype=true)) + 
    geom_line() + 
    geom_vline(xintercept = c(mod$lambda.min, mod$lambda.1se),linetype="dashed")+
    annotate(geom = "text",
             label = c("lambda.min","lambda.1se",mod_df_min, mod_df_1se),
             x = c(mod$lambda.min, mod$lambda.1se,mod$lambda.min, mod$lambda.1se),
             y = c(-0.10,-0.10, 0.2,0.2),
             angle = 90, 
             vjust = 1,
             size=3)+
    scale_x_log10() + 
    ylim(-0.20, 0.25)+
    xlab("Lambda (log scale)") + 
    guides(color = guide_legend(title = ""), 
           linetype = guide_legend(title = "")) +
    ggtitle(mod_title)+
    theme_bw() + 
    theme(legend.key.width = unit(3,"lines"))
}




## generate gglasso model under different standardization for plot
plot_data <- function(data_rda){
  
  # load X, Z
  load(data_rda)
  
  # true coefficients
  beta_0 <- 0
  beta_p <- c(rep(0.05,3),rep(-0.05,3),rep(0,6))
  
  # true model
  linear_part <- beta_0 + (as.matrix(X) %*% beta_p)
  Y <- rbinom(nrow(X),size = 1, prob = ilogit(linear_part))
  Y_np1 <- ifelse(Y==0,-1,1)
  
  # group index
  group_z <- c(c(1:4),rep(5,2),rep(6,3),c(7:10),rep(11,2),rep(12,3))
  
  # 1. none standardization
  dummies1 <- dummyVars(~., data = Z, fullRank = T)
  Z_dum1 <- predict(dummies1,Z)
  
  # 2. Andrew: scale continuous by 2sd
  Z_dum2 <- Z_dum1 %>% as.data.frame() %>% 
    mutate_at(c('Z1','Z7'),function(x) return(x/2/sd(x))) %>% 
    as.matrix()
  
  # 3. Min-max continuous
  Z_dum3 <- Z_dum1 %>% as.data.frame() %>% 
    mutate_at(c('Z1','Z7'),function(x) return((x-min(x))/(max(x)-min(x)))) %>% 
    as.matrix()
  
  # 4. dichotomize continuous, and reduce multicatgory
  Z_bi4 <- Z %>% 
    mutate(Z1=ifelse(Z1<=quantile(Z1,probs = 0.5),0,1)) %>% 
    mutate(Z7=ifelse(Z7<=quantile(Z7,probs = 0.5),0,1)) %>% 
    mutate_at(c('Z5','Z6'), function(x) multi_fix(x)) %>%
    mutate_at(c('Z11','Z12'), function(x) multi_random(x)) %>%as.matrix()
  
  # 5. Proposed: dichotomize continuous, and reduce multicatgory, standard all binary by 2sd
  Z_bi5 <- Z_bi4 %>% as.data.frame() %>%  
    mutate_all(., function(x) return(x/2/sd(x))) %>% 
    as.matrix()
  
  # 6. Z-score: stand all by 1sd
  Z_dum6 <- Z_dum1 %>% as.data.frame() %>% 
    mutate_all(.,function(x) return(x/sd(x))) %>% 
    as.matrix()
  
  # cv folder number
  K <- 5
  
  mod_1 <- cv.gglasso(x=Z_dum1,y=Y_np1,group=group_z,loss="logit",pred.loss = "misclass",nfolds = K)
  
  mod_2 <- cv.gglasso(x=Z_dum2,y=Y_np1,group=group_z,loss="logit",pred.loss = "misclass",nfolds = K)
  
  mod_3 <- cv.gglasso(x=Z_dum3,y=Y_np1,group=group_z,loss="logit",pred.loss = "misclass",nfolds = K)
  
  mod_4 <- cv.gglasso(x=Z_bi4,y=Y_np1,group=c(1:ncol(X)),loss="logit",pred.loss = "misclass",nfolds = K)
  
  mod_5 <- cv.gglasso(x=Z_bi5,y=Y_np1,group=c(1:ncol(X)),loss="logit",pred.loss = "misclass",nfolds = K)
  
  mod_6 <- cv.gglasso(x=Z_dum6,y=Y_np1,group=group_z,loss="logit",pred.loss = "misclass",nfolds = K)
  
  return(list(mod_1=mod_1,mod_2=mod_2,mod_3=mod_3,mod_4=mod_4,mod_5=mod_5,mod_6=mod_6))
}




# Transform X to Z: normal distribution
Z_gen <- function(X){

  Z <- data.frame(Z1=X$X1,
                  Z2=ifelse(X$X2<=qnorm(p=0.5, mean=0, sd=1),0,1),
                  Z3=ifelse(X$X3<=qnorm(p=0.7, mean=0, sd=1),0,1),
                  Z4=ifelse(X$X4<=qnorm(p=0.9, mean=0, sd=1),0,1),
                  Z5=factor(ifelse(X$X5<=qnorm(p=1/3, mean=0, sd=1),"f0",
                                   ifelse(X$X5<=qnorm(p=2/3, mean=0, sd=1),"f1","f2")),
                            levels = paste0("f",c(0,1,2))),
                  Z6=factor(ifelse(X$X6<=qnorm(p=1/4, mean=0, sd=1),"f0",
                                   ifelse(X$X6<=qnorm(p=2/4, mean=0, sd=1),"f1",
                                          ifelse(X$X6<=qnorm(p=3/4, mean=0, sd=1),"f2","f3"))),
                            levels = paste0("f",c(0,1,2,3))),
                  Z7=X$X7,
                  Z8=ifelse(X$X8<=qnorm(p=0.5, mean=0, sd=1),0,1),
                  Z9=ifelse(X$X9<=qnorm(p=0.7, mean=0, sd=1),0,1),
                  Z10=ifelse(X$X10<=qnorm(p=0.9, mean=0, sd=1),0,1),
                  Z11=factor(ifelse(X$X11<=qnorm(p=1/3, mean=0, sd=1),"f0",
                                    ifelse(X$X11<=qnorm(p=2/3, mean=0, sd=1),"f1","f2")),
                             levels = paste0("f",c(0,1,2))),
                  Z12=factor(ifelse(X$X12<=qnorm(p=1/4, mean=0, sd=1),"f0",
                                    ifelse(X$X12<=qnorm(p=2/4, mean=0, sd=1),"f1",
                                           ifelse(X$X12<=qnorm(p=3/4, mean=0, sd=1),"f2","f3"))),
                             levels = paste0("f",c(0,1,2,3)))

  )

  return(Z)
}




# Transform X to Z: exp distribution
Z_gen_exp <- function(X){
  Z <- data.frame(Z1=X$X1,
                  Z2=ifelse(X$X2<=qexp(p=0.5,rate=1),0,1),
                  Z3=ifelse(X$X3<=qexp(p=0.7,rate=1),0,1),
                  Z4=ifelse(X$X4<=qexp(p=0.9,rate=1),0,1),
                  Z5=factor(ifelse(X$X5<=qexp(p=1/3,rate=1),"f0",
                                   ifelse(X$X5<=qexp(p=2/3,rate=1),"f1","f2")),
                            levels = paste0("f",c(0,1,2))),
                  Z6=factor(ifelse(X$X6<=qexp(p=1/4,rate=1),"f0",
                                   ifelse(X$X6<=qexp(p=2/4,rate=1),"f1",
                                          ifelse(X$X6<=qexp(p=3/4,rate=1),"f2","f3"))),
                            levels = paste0("f",c(0,1,2,3))),
                  Z7=X$X7,
                  Z8=ifelse(X$X8<=qexp(p=0.5,rate=1),0,1),
                  Z9=ifelse(X$X9<=qexp(p=0.7,rate=1),0,1),
                  Z10=ifelse(X$X10<=qexp(p=0.9,rate=1),0,1),
                  Z11=factor(ifelse(X$X11<=qexp(p=1/3,rate=1),"f0",
                                    ifelse(X$X11<=qexp(p=2/3,rate=1),"f1","f2")),
                             levels = paste0("f",c(0,1,2))),
                  Z12=factor(ifelse(X$X12<=qexp(p=1/4,rate=1),"f0",
                                    ifelse(X$X12<=qexp(p=2/4,rate=1),"f1",
                                           ifelse(X$X12<=qexp(p=3/4,rate=1),"f2","f3"))),
                             levels = paste0("f",c(0,1,2,3)))
                  
  )
  return(Z)
}




# Transform X to Z: gamma distribution
Z_gen_gamma <- function(X){
  
  g_alpha <- 10
  g_beta <- 0.32
  
  Z <- data.frame(Z1=X$X1,
                  Z2=ifelse(X$X2<=qgamma(p = 0.5, shape =g_alpha, scale = g_beta),0,1),
                  Z3=ifelse(X$X3<=qgamma(p = 0.7, shape =g_alpha, scale = g_beta),0,1),
                  Z4=ifelse(X$X4<=qgamma(p = 0.9, shape =g_alpha, scale = g_beta),0,1),
                  Z5=factor(ifelse(X$X5<=qgamma(p = 1/3, shape =g_alpha, scale = g_beta),"f0",
                                   ifelse(X$X5<=qgamma(p = 2/3, shape =g_alpha, scale = g_beta),"f1","f2")),
                            levels = paste0("f",c(0,1,2))),
                  Z6=factor(ifelse(X$X6<=qgamma(p = 1/4, shape =g_alpha, scale = g_beta),"f0",
                                   ifelse(X$X6<=qgamma(p = 2/4, shape =g_alpha, scale = g_beta),"f1",
                                          ifelse(X$X6<=qgamma(p = 3/4, shape =g_alpha, scale = g_beta),"f2","f3"))),
                            levels = paste0("f",c(0,1,2,3))),
                  Z7=X$X7,
                  Z8=ifelse(X$X8<=qgamma(p = 0.5, shape =g_alpha, scale = g_beta),0,1),
                  Z9=ifelse(X$X9<=qgamma(p = 0.7, shape =g_alpha, scale = g_beta),0,1),
                  Z10=ifelse(X$X10<=qgamma(p = 0.9, shape =g_alpha, scale = g_beta),0,1),
                  Z11=factor(ifelse(X$X11<=qgamma(p = 1/3, shape =g_alpha, scale = g_beta),"f0",
                                    ifelse(X$X11<=qgamma(p = 2/3, shape =g_alpha, scale = g_beta),"f1","f2")),
                             levels = paste0("f",c(0,1,2))),
                  Z12=factor(ifelse(X$X12<=qgamma(p = 1/4, shape =g_alpha, scale = g_beta),"f0",
                                    ifelse(X$X12<=qgamma(p = 2/4, shape =g_alpha, scale = g_beta),"f1",
                                           ifelse(X$X12<=qgamma(p = 3/4, shape =g_alpha, scale = g_beta),"f2","f3"))),
                             levels = paste0("f",c(0,1,2,3)))
                  
  )
  
  return(Z)
}
