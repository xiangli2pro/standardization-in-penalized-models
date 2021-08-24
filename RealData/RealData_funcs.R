##---------------------------------------------------------------

## Functions called in RealData_models.R and RealData_models_eval.R

## Author: Xiang Li
## Date: 03-17-2021

## Functions created:
## LRT_func(), multi_reduce(), coef.rep(), cv_wgglasso()
## pattern_wgglasso()
## se_bic_func(), se_auc_func(), get_features(), cbind.fill(), get_form()

##---------------------------------------------------------------




##----------------
## LRT_func
##----------------
## likelihood-ratio test of multi-category variables
LRT_func <- function(x,y){
  mod <- data.frame(y,x) %>% glm(y~x, data = .,family = "binomial") 
  test <- anova(mod,test = "LRT")
  test[2,5]<0.05
}




##----------------
## multi_reduce
##----------------
## function: 
## dichotomize continuous variable to binary
## regroup the multi-category to binary
multi_reduce <- function(x,y=Y, p=0.5, q=0.5){
  
  # parameters
  # p: cutpoints of continuous covariate
  # q: cutpoints of multi-category covariate
  
  # dichotomize continuous by p-th quantile
  if (length(levels(x))==0){
    
    new_x <- ifelse(x<=quantile(x,probs = p),0,1)
    return(new_x)
    
  # turn binary to numerical 0, 1
  } else if (length(levels(x))==2){
    
    g0_index <- (x == levels(x)[1])
    new_x <- vector(length = length(x))
    new_x[g0_index] <- 0;
    new_x[!g0_index] <- 1;
    
    return(new_x)
    
  # regroup the multi-category
  } else if (length(levels(x))>2){
    
    x_levels <- levels(x)
    
    # LRT test
    if(LRT_func(x,y)){
      
      mod_coef <- data.frame(y,x) %>% 
        glm(y~x, data = .,family = "binomial") %>% coef()
      mod_coef[2:length(mod_coef)] <- mod_coef[2:length(mod_coef)] + mod_coef[1] 
      
      # rank coefficients in ascending order
      x_coef_order <- order(mod_coef)  
      
    } else {
      
      # if LRT-test is not significant, randomly group the multi-category
      x_coef_order <- sample(1:length(x_levels),length(x_levels),replace = FALSE)
    }
    
    # frequency of each level by the order of coefficient size
    x_freq_order <- table(x)[x_coef_order]/length(x) 
    
    # find the  level with corresponding cumulative frequency greater than q
    potential_cut <- which(cumsum(x_freq_order)>=q)[1] 
    
    # adjust the potential_cut level to make the binary prob. closest to q
    # e.g levels: c(1, 2, 3, 4) with cumulative freq. c(0.1, 0.45, 0.6, 1), q=0.5
    # we should set potential_cut = 2 not 3, because 0.45 is closer to 0.5 than 0.6
    if (potential_cut>1){
      
      # avoid the case of choosing level with cumulative freq. 1 as cutoff level
      if(potential_cut==length(x_levels)){
        
        x_cut <- potential_cut-1
    
      } else {
        
        # test which level have cumulative frequency closer to q
        lower_diff <- abs(q-cumsum(x_freq_order)[potential_cut-1])
        upper_diff <- abs(q-cumsum(x_freq_order)[potential_cut])
        x_cut <- ifelse(lower_diff<upper_diff, potential_cut-1, potential_cut)
      }
      
    } else {
      
      x_cut <- potential_cut
    }
    
    g0_index <- x %in% (x_levels[x_coef_order])[c(1:x_cut)]
    
    new_x <- vector(length = length(x))
    new_x[g0_index] <- 0 
    new_x[!g0_index] <- 1 
    
    return(new_x)
  } 
}




##----------------
## coef.rep
##----------------
## get group index for all the dummy covariates
coef.rep <- function(X.dt){
  
  # X.dt is the original dataset
  coef.rep.ind <- c()
  
  for( i in 1: ncol(X.dt)){
    
    if ( is.factor(X.dt[,i]) == TRUE){
      rep.num <- length(levels(X.dt[,i]))-1 
      coef.rep.ind <- c(coef.rep.ind,rep(i,rep.num))
      
    }else{
      
      coef.rep.ind<-c(coef.rep.ind,i)
    }
  }
  
  return(coef.rep.ind)
}




##----------------
## cv_wgglasso
##----------------
## cross-validation performance of gglasso on data standardized by different methods
cv_wgglasso <- function(y, X_dum, X_weight, K,K_folds, coef_group, grp_lambda_seq){
  
  # train and test data set
  train_index <- do.call(c,K_folds[-K])
  test_index <- unlist(K_folds[K])
  
  train_data <- list(y=y[train_index],
                     X_weight=X_weight[train_index],
                     X=cbind(1,X_dum[train_index,]))
  
  test_data <- list(y=y[test_index],
                    X_weight=X_weight[test_index],
                    X=cbind(1,X_dum[test_index,]))
  
  # fit grplasso on train data
  train_fit <- grplasso(x=train_data$X, y= train_data$y, index=coef_group, weights = train_data$X_weight, 
                        lambda = grp_lambda_seq,
                        model = LogReg(), penscale = sqrt, 
                        center = FALSE, standardize = FALSE)
  
  # predict test data
  test_fit <- predict(train_fit, newdata=test_data$X, type="response")
  
  # extract AUC score of test prediction
  auc_func <- function(x){
    auc(roc(test_data$y,x))
  }
  
  auc_score <- apply(test_fit,2,auc_func)

  return(auc_score)
}




##----------------
## pattern_wgglasso
##----------------
## evaluate model performance of grplasso on different standardized data
pattern_wgglasso <- function(X, y, X_weight, grp_lambda_seq, p=0.9, q=0.9, stand){
  
  # identity covariates names of different functional forms
  X_type <- unlist(lapply(X, function(x) length(levels(x))))
  X_cont <- names(X)[which(X_type==0)]
  X_bi <- names(X)[which(X_type==2)]
  X_multi <- names(X)[which(X_type>2)]
  
  # standardize data
  if (stand == "none"){
    
    dummies <- dummyVars(~., data = X, fullRank = T)
    X_dum <- predict(dummies,X) 
    
    coef_group <- c(NA,coef.rep(X))
    
  } else if (stand == "zscore"){
    
    dummies <- dummyVars(~., data = X, fullRank = T)
    X_dum <- scale(predict(dummies,X)) 
    
    coef_group <- c(NA,coef.rep(X))
    
  } else if (stand == "gelman"){
    
    scale2_func <- function(x,na.rm=TRUE){
      x/2/sd(x,na.rm)
    }
    
    X_2sd <- X %>% mutate_at(c("AGE","numother","TIMEMD"),scale2_func)
    dummies <- dummyVars(~., data = X_2sd, fullRank = T)
    X_dum <- predict(dummies,X_2sd)
    
    coef_group <- c(NA,coef.rep(X))
    
  } else if (stand == "minmax"){
    
    min_max_func <- function(x){
      (x-min(x))/(max(x)-min(x))
    }
    
    X_mm <- X %>% mutate_at(c("AGE","numother","TIMEMD"),min_max_func)
    dummies <- dummyVars(~., data = X_mm, fullRank = T)
    X_dum <- predict(dummies,X_mm)
    
    coef_group <- c(NA,coef.rep(X))
    
  } else if (stand == "propose"){
    
    X_dum <- lapply(X, function(x) multi_reduce(x,y=y,p=p,q=q)) %>%
      as.data.frame() %>%
      mutate_all(., function(x) x/2/sd(x)) %>% 
      as.matrix()
    
    coef_group <- c(NA,c(1:ncol(X)))
  }
  
  # fit model
  mod_fit <- grplasso(x=cbind(1,X_dum), y= y, index=coef_group, weights =X_weight, 
                      lambda = grp_lambda_seq,
                      model = LogReg(), penscale = sqrt, 
                      center = FALSE, standardize = FALSE)
  
  # get coefficients
  coef_fit <- mod_fit$coefficients
  coef_idx <- alply(coef_fit,2, function(x) which(x!=0)) # use alply to return list
  coef_names <- lapply(coef_idx, function(x) (names(X)[unique(coef_group[x])])[-1])
  
  # refit logistic regression using the selected covariates
  # get AIC, BIC, log-likelihood, AUC of each model
  # get no. of selected continuous, selected binary, selected mutli-category
  eval_func <- function(names_select, y){
    
    if (length(names_select)==0){
      
      return(rep(NA,7))
      
    } else {
      
      mod_glm <- glm(y~., data=data.frame(y,X[,names_select]), family = "binomial", weights = X_weight)
      return(c(logLik(mod_glm), AIC(mod_glm), BIC(mod_glm), auc(roc(y,mod_glm$fitted.values)), 
               sum(names_select %in% X_cont), sum(names_select %in% X_bi), sum(names_select %in% X_multi)))
    }
  }
  
  eval_mat <- vapply(coef_names,function(x) eval_func(x,y), numeric(7))
  
  return(eval_mat)
}




##----------------
## se_bic_func
##----------------
## find bic_1se
se_bic_func <- function(x){
  min_val <- min(x,na.rm = TRUE)
  se_val <- sd(x,na.rm = TRUE)
  return(which.min(abs(x-(min_val+se_val))))
}




##----------------
## se_auc_func
##----------------
## find auc_1se
se_auc_func <- function(x){
  max_val <- max(x,na.rm = TRUE)
  se_val <- sd(x,na.rm = TRUE)
  return(which.min(abs(x-(max_val-se_val))))
}




##----------------
## cbind.fill
##----------------
## 
cbind.fill <- function(...){
  nm <- list(...) 
  nm <- lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}




##----------------
## get_form
##----------------
## get functional form of a covariate
get_form <- function(x){
  if (x %in% X_bi){
    return('binary')
  } else if (x %in% X_multi){
    return('multi')
  }
}





##----------------
## get_features
##----------------
## return names of selected covariates
get_features <- function(X, y, X_weight, lambda, p=0.9, q=0.9, stand){
  
  X_type <- unlist(lapply(X, function(x) length(levels(x))))
  X_cont <- names(X)[which(X_type==0)]
  X_bi <- names(X)[which(X_type==2)]
  X_multi <- names(X)[which(X_type>2)]
  
  if (stand == "none"){
    
    dummies <- dummyVars(~., data = X, fullRank = T)
    X_dum <- predict(dummies,X) ## none
    
    coef_group <- c(NA,coef.rep(X))
    
  } else if (stand == "zscore"){
    
    dummies <- dummyVars(~., data = X, fullRank = T)
    X_dum <- scale(predict(dummies,X)) ## z-score
    
    coef_group <- c(NA,coef.rep(X))
    
  } else if (stand == "gelman"){
    
    scale2_func <- function(x,na.rm=TRUE){
      x/2/sd(x,na.rm)
    }
    
    X_2sd <- X %>% mutate_at(c("AGE","numother","TIMEMD"),scale2_func)
    dummies <- dummyVars(~., data = X_2sd, fullRank = T)
    X_dum <- predict(dummies,X_2sd)
    
    coef_group <- c(NA,coef.rep(X))
    
  } else if (stand == "minmax"){
    
    min_max_func <- function(x){
      (x-min(x))/(max(x)-min(x))
    }
    
    X_mm <- X %>% mutate_at(c("AGE","numother","TIMEMD"),min_max_func)
    dummies <- dummyVars(~., data = X_mm, fullRank = T)
    X_dum <- predict(dummies,X_mm)
    
    coef_group <- c(NA,coef.rep(X))
    
  } else if (stand == "propose"){
    
    X_dum <- lapply(X, function(x) multi_reduce(x,y=y,p=p,q=q)) %>% as.data.frame() %>%
      mutate_all(., function(x) x/2/sd(x)) %>% 
      as.matrix()
    
    coef_group <- c(NA,c(1:ncol(X)))
  }
  
  mod_fit <- grplasso(x=cbind(1,X_dum), y= y, index=coef_group, weights =X_weight, 
                      lambda = lambda,
                      model = LogReg(), penscale = sqrt, 
                      center = FALSE, standardize = FALSE)
  
  coef_fit <- mod_fit$coefficients
  coef_idx <- which(coef_fit!=0)
  coef_names <- (names(X)[unique(coef_group[coef_idx])])[-1]
  
  return(coef_names)
}


