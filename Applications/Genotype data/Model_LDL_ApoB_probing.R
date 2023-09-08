################################################################
################################################################
# ---------- Application genotype data LDL and ApoB  --------- #

library(gamboostLSS)
library(gamlss.dist)
library(copula)

source("Copulas/Copula_Gaussian.R")
source("Copulas/Copula_Clayton.R")
source("Copulas/Copula_Gumbel.R")

sim <- function(seed){
  
  load("data_LDL_adj_ApoB.RData")
  data_chol <- data[[1]]$LDL_adj
  
  load("data_ApoB_adj_LDL.RData")
  dataset <- data[[1]]
  weights <- data[[2]]
  dataset <- cbind(LDL_adj = data_chol,dataset)
  train = dataset

  weights_train=weights
  train <- train[weights_train == 1,]

  set.seed(seed)

  probes.train <- as.data.frame(sapply(train[,-c(1,2)], sample)) # shuffling all variables (leaving out response which is at position 1 and 2)
  names(probes.train) <- paste(names(train)[-c(1,2)], "probe", sep = "_")
  probes.train <- cbind(train,probes.train)
  
  load("data_test_LDL_adj_ApoB.RData")
  test = test 
  probes.test <- data.frame(matrix(0,ncol = dim(test)[2],nrow = n))
  names(probes.test) <- paste(names(test)[-c(1,2)], "probe", sep = "_")

  time_classic = Sys.time()
  mod = glmboostLSS(cbind(LDL_adj,ApoB) ~., data = probes.train, control = boost_control(mstop = 1000, nu = 0.01,trace = T), method = 'noncyclic', families =  Gauss_Cop(marg1 = "LOGLOG", marg2 = "GA", stabilization = "MAD"))
 
  
  p  = ncol(train)- 2
  if(!any(selected(mod,merge = T) > p)) mod[3000]
  if(!any(selected(mod,merge = T) > p)) mod[5000]
  mstop_probes <- min(which(selected(mod,merge = T) > p   )) - 1 
  mod[mstop_probes]
  runtime_classic = Sys.time() - time_classic
  units(runtime_classic) = "mins"
  
  mstop.mod<-  vector('list')
  mstop.mod$mstop <- mstop_probes
  mstop.mod$mu1 <- mod$mu1$mstop()
  mstop.mod$mu2 <- mod$mu2$mstop()
  mstop.mod$sigma1 <- mod$sigma1$mstop()
  mstop.mod$sigma2 <- mod$sigma2$mstop()
  mstop.mod$rho <- mod$rho$mstop()
  
  Coefficients <- coef(mod)
  
  probes.test <- cbind(test, probes.test)
  #Prediction error
  pred_mu1    <- predict(mod$mu1, newdata = probes.test, type = 'response')
  pred_sigma1 <- predict(mod$sigma1, newdata = probes.test, type = 'response')
  pred_mu2    <- predict(mod$mu2, newdata = probes.test, type = 'response')
  pred_sigma2 <- predict(mod$sigma2, newdata = probes.test, type = 'response')
  pred_rho    <- predict(mod$rho, newdata = probes.test, type = 'response')
  
  pred <- vector('list')
  pred$mu1 <- pred_mu1
  pred$mu2 <- pred_mu2
  pred$sigma1 <- pred_sigma1
  pred$sigma2 <- pred_sigma2
  pred$rho <- pred_rho
  
  ### helpers LOGLOG
  ### pdf 
  dLogLogistic <- function(x, mu = 1, sigma = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(sigma < 0)) 
      stop(paste("sigma must be positive", "\n", ""))
    
    
    z <- (x/mu)^(sigma)
    
    loglik <- log(sigma) + (sigma) * (log(x) - log(mu)) - log(x) - 2 * log(1 + z)
    
    if (log == FALSE) 
      ft <- exp(loglik)
    else ft <- loglik
    ft
  }
  
  ### cdf 
  
  pLogLogistic <- function(x, mu = 1, sigma = 1, log = FALSE){
    
    if (any(mu < 0)) 
      stop(paste("mu must be positive", "\n", ""))
    
    if (any(sigma < 0)) 
      stop(paste("sigma must be positive", "\n", ""))
    
    
    log_cdf <- log(1/(1 + (x/mu)^(-sigma)))
    
    if (log == FALSE) 
      ft <- exp(log_cdf)
    else ft <- log_cdf
    ft
  }
  
  
  
  loss <- function(y, mu1, sigma1, mu2, sigma2, rho){-(log(dLogLogistic(x = y[, 1], mu = mu1, sigma = sigma1)) +  
                                                         log(dgamma(x = y[, 2], shape = 1/sigma2^2, scale = mu2 *  
                                                                      sigma2^2)) - 0.5 * log(1 - rho^2) + (rho/(1 - rho^2)) *  
                                                         qnorm(pLogLogistic(x = y[, 1], mu = mu1, sigma = sigma1)) *  
                                                         qnorm(pgamma(q = y[, 2], shape = 1/sigma2^2, scale = mu2 *  
                                                                        sigma2^2)) - (rho^2/(2 * (1 - rho^2))) * (qnorm(pLogLogistic(x = y[,  
                                                                                                                                           1], mu = mu1, sigma = sigma1))^2 + qnorm(pgamma(q = y[,  
                                                                                                                                                                                                 2], shape = 1/sigma2^2, scale = mu2 * sigma2^2))^2)) 
  }
  
  
  LogLik_values <- loss(y = probes.test[,c('LDL_adj','ApoB')], mu1 = pred_mu1, mu2 = pred_mu2,
                        sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho)
  LogLik <- sum(loss(y = probes.test[,c('LDL_adj','ApoB')], mu1 = pred_mu1, mu2 = pred_mu2,
                     sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho),na.rm=T)


  return(list(Coefficients=Coefficients, mstop=mstop.mod, pred=pred, LogLik_values=LogLik_values, LogLik=LogLik))

}

results = mclapply(1, sim, mc.cores = 1, mc.preschedule = FALSE)

save(results, file="LDL_ApoB_models_prob.RData")
