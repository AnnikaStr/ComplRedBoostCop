###########################################################
###########################################################
# ------ Simulation Gumbel - Modification of Hans et al-- #

library("BatchJobs")
library("Rcpp")
library("gamboostLSS")
library("gamlss.dist")
library("copula")
library("scoringRules")

source("defs.R")
source("Copulas/Copula_Gaussian.R")
source("Copulas/Copula_Clayton.R")
source("Copulas/Copula_Gumbel.R")
source("DeselectBoost_Update.R")

reg <- makeExperimentRegistry("VariableSelection_Gumbel", seed = 123L, packages = c("gamboostLSS","copula"))
reg$default.resources = list(memory = MEMORY, ntasks = NTASKS, walltime = WALLTIME)
reg$max.concurrent.jobs = 1000L
reg$cluster.functions = makeClusterFunctionsMulticore(ncpus = NCPUS)

dat <- function(data,job){
  
  p = data$p
  n_train = data$n
  n_test  = data$n_test
  n_eval = data$n_eval
  
  
  ### Data generating process
  # non-linear effects
  s1 <- function(x) 0.5 * cos(pi * x)
  s2 <- function(x) 1 - exp((x + 1) * 2) * 0.02
  s3 <- function(x) 0.5 * sin(pi * x)
  s4 <- function(x) 1.5 * log(4.5 - 1.7 * sin(pi * x))
  
  # Simulation function (see: A Bivariate Copula Additive Model for Location, Scale and Shape (2017), Marra and Radice)
  data.gen <- function(FAM, s1, s2, s3, s4, Z, i){
    
    
    z1 <- Z[i, 1]
    z2 <- Z[i, 2]
    z3 <- Z[i, 3]
    z4 <- Z[i, 4]
    
    eta_mu1 <- -0.75 * z1 + s1(z3)
    eta_sigma1 <- -0.7
    eta_mu2 <- - 0.5 - 0.7 * z1 + s2(z2)
    eta_sigma2 <- 2 + 0.5 * z2
    eta_theta <- -0.8 + s4(z4)
    
    if( FAM == "clayton") theta.para <- exp(eta_theta) + 1e-07
    if( FAM == "gumbel") theta.para <- exp(eta_theta) + 1
    if( FAM %in% c ("normal" , "t")) theta.para <- tanh(eta_theta)              # This might be the problem! other link function!!!!
    if( FAM %in% c ("clayton" , "gumbel") ) Cop <- archmCopula(family = FAM,
                                                               dim = 2, param = theta.para)
    else Cop <- ellipCopula(family = FAM, dim = 2, param = theta.para , df = 4)
    speclist1 <- list( meanlog = eta_mu1, sdlog = exp(eta_sigma1) )
    speclist2 <- list( mu = exp(eta_mu2), sigma = exp(eta_sigma2), nu = 1, tau = 1)
    spec <- mvdc(copula = Cop, c("lnorm" , "GB2") , list(speclist1 , speclist2) )
    c(rMvdc(1, spec), Z[i, ])
  }
  
  n_eval = 1500
  
  data_train <- as.data.frame(matrix(0, ncol = p+2+1, nrow=(n_train+n_eval)))
  names(data_train) <- c("y1", "y2", paste0("x", 1:p),'Intercept')

  X_train <- matrix(NA, ncol = p, nrow = (n_train+n_eval))
  
  for(i in 1:p){
    X_train[,i] <- runif((n_train+n_eval), -1, 1)
  }
  
  X_train <- cbind(X_train, 1)
  
  weights_train <- c(rep(1, times = n_train), rep(0, times = n_eval))
  
  for(i in 1:(n_train + n_eval)){
    data_run <- data.gen(FAM = "gumbel", s1 = s1, s2 = s2, s3 = s3, s4 = s4, Z = X_train, i = i)
    data_train[i, ] <- data_run
  }
  
  n_test = 1000
  ### generate test data
  data_test <- as.data.frame(matrix(0, ncol = 2 + p +1, nrow = (n_test)))
  names(data_test) <- c("y1", "y2", paste0("x", 1:p),'Intercept')
  X_test <- (matrix(0, ncol = p, nrow = n_test))
  
  for(i in 1:p){
    X_test[,i] <- runif((n_test), -0.95, 0.95)    
  }
  
  X_test <- cbind(X_test, 1)
  
  for(i in 1:n_test){
    data_run <- data.gen(FAM = "gumbel", s1 = s1, s2 = s2, s3 = s3, s4 = s4, Z = X_test, i = i)
    data_test[i, ] <- data_run
  }
  
  TrueBeta <-  vector('list')
  TrueBeta$mu1 <- c(1,0,1,rep(0,times = p-3))
  TrueBeta$sigma1 <- c(0,0,0,rep(0,times = p-3))
  TrueBeta$mu2 <- c(1,1,0,rep(0,times = p-3))
  TrueBeta$sigma2 <- c(0,1,0,rep(0,times = p-3))
  TrueBeta$rho <- c(0,0,0,1,rep(0,times = p-4))
  
  seed = runif(1,1,100)
  
  list(TrueBeta = TrueBeta, train = data_train, test = data_test, weights_train=weights_train, p = p, n = n_train, seed = seed)
  
}


classic <- function(job,data,instance){

  train <- instance$train
  test <- instance$test
  TrueBeta <- instance$TrueBeta
  weights_train <- instance$weights_train
  seed = instance$seed
  
  set.seed(seed)
  
  form <- as.formula(( paste( "cbind(y1,y2) ~",  paste(c(paste("bbs(x", 1:p, ")", sep = "")), collapse = "+"))))
  
  time = Sys.time()
  mod = gamboostLSS(form, data = train, control = boost_control(mstop = 5000, nu = 0.01, risk = "oobag"), weights = weights_train, method = 'noncyclic', families = Gumbel_Cop(marg1 = "LOGNO",
                                                                                                                                                                               marg2 = "LOGLOG"))
  MSTOP <- which.min(risk(mod,merge = T))
  oobag.risk <- risk(mod,merge = T)
  
  rm(mod) # removed the fist fitted model
  train <- train[weights_train == 1, ]
  
  mod = gamboostLSS(form, data = train, control = boost_control(mstop = MSTOP, nu = 0.01), method = 'noncyclic', families = Gumbel_Cop(marg1 = "LOGNO",marg2 = "LOGLOG"))
  runtime = Sys.time() - time
  units(runtime) = "mins"
  
  mstop.mod<-  vector('list')
  mstop.mod$mstop <- MSTOP
  mstop.mod$mu1 <- mod$mu1$mstop()
  mstop.mod$mu2 <- mod$mu2$mstop()
  mstop.mod$sigma1 <- mod$sigma1$mstop()
  mstop.mod$sigma2 <- mod$sigma2$mstop()
  mstop.mod$rho <- mod$rho$mstop()
  
  Coefficients <- coef(mod)
  
  TP <- vector(length = 5, 'list')
  TP <- lapply(1:5,function(i){sum( gsub("^[^(]*\\(([^)]*)\\)","\\1", names(Coefficients[[i]])) %in% paste('x',which(TrueBeta[[i]]!=0),sep = ''))})
  names(TP) <- names(Coefficients)
  
  FP <- vector(length = 5, 'list')
  FP <- lapply(1:5,function(i){sum(!( gsub("^[^(]*\\(([^)]*)\\)","\\1", names(Coefficients[[i]])) %in% paste('x',which(TrueBeta[[i]]!=0),sep = '')))})
  names(FP) <- names(Coefficients)
  
  #Prediction error
  pred_mu1    <- predict(mod$mu1, newdata = test, type = 'response')
  pred_sigma1 <- predict(mod$sigma1, newdata = test, type = 'response')
  pred_mu2    <- predict(mod$mu2, newdata = test, type = 'response')
  pred_sigma2 <- predict(mod$sigma2, newdata = test, type = 'response')
  pred_rho    <- predict(mod$rho, newdata = test, type = 'response')
  
  pred <- vector('list')
  pred$mu1 <- pred_mu1
  pred$mu2 <- pred_mu2
  pred$sigma1 <- pred_sigma1
  pred$sigma2 <- pred_sigma2
  pred$rho <- pred_rho
  
  ### helpers LOGLOG
  ### pdf 
  dLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
    
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
  pLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
    
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
  
  loss <- function (y, mu1, sigma1, mu2, sigma2, rho){
    -(log(dlnorm(x = y[, 1], meanlog = mu1, sdlog = sigma1)) +
        log(dLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2)) +
        log(rho + 1) - (1 + rho) * (log(plnorm(q = y[, 1], meanlog = mu1, sdlog = sigma1)) +
                                      log(pLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2))) -
        (2 + 1/rho) * log(plnorm(q = y[, 1], meanlog = mu1, sdlog = sigma1)^(-rho) +
                            pLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2)^(-rho) - 1))}
  
  LogLik_values <- loss(y = test[,c('y1','y2')], mu1 = pred_mu1, mu2 = pred_mu2,
                        sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho)
  LogLik <- sum(loss(y = test[,c('y1','y2')], mu1 = pred_mu1, mu2 = pred_mu2,
                     sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho),na.rm=T)
  
  es_copula <- rep(NA, 1000)
  
  for(pred_sim_run in 1:1000){
    
    pred_sample_copula <- matrix(NA, nrow = 2, ncol = 10000)
    
    Cop <- archmCopula(family = "gumbel", dim = 2, param = pred_rho[pred_sim_run])
    speclist1 <- list( meanlog = pred_mu1[pred_sim_run], sdlog = pred_sigma1[pred_sim_run] )
    speclist2 <- list( mu = pred_mu2[pred_sim_run], sigma = pred_sigma2[pred_sim_run], nu = 1, tau = 1)
    spec <- mvdc(copula = Cop, c("lnorm" , "GB2") , list(speclist1 , speclist2) )
    samp_cop <-  rMvdc(10000, spec)
    pred_sample_copula[1, ] <- samp_cop[ ,1]
    pred_sample_copula[2, ] <- samp_cop[ ,2]
    
    es_copula[pred_sim_run] <- es_sample(y = c(test$y1[pred_sim_run], test$y2[pred_sim_run]), dat = pred_sample_copula)
    
  }
  list(instance, Coefficients=Coefficients, TP=TP, FP=FP, mstop=mstop.mod, pred=pred, LogLik_values=LogLik_values, LogLik=LogLik,  energy_score =es_copula, runtime=runtime)
}



probing <- function(job, data, instance){

  train <- instance$train
  weights_train <- instance$weights_train
  train <- train[weights_train == 1,]
  test <- instance$test
  TrueBeta <- instance$TrueBeta
  seed = instance$seed
  
  set.seed(seed)
  
  probes.train <- as.data.frame(sapply(train[,-c(1,2,p+2+1)], sample)) # shuffling all variables (leaving out response which is at position 1 and 2)
  names(probes.train) <- paste(names(train)[-c(1,2,p+2+1)], "probe", sep = "_")
  probes.train <- cbind(train,probes.train)
  
  probes.test <- data.frame(matrix(0,ncol = p+1,nrow = n))
  names(probes.test) <- paste(names(test)[-c(1,2)], "probe", sep = "_")
  probes.test <- cbind(test, probes.test)

  time = Sys.time()
  mod = gamboostLSS(cbind(y1,y2) ~., data = probes.train, control = boost_control(mstop = 2000, nu=0.01),
                    method = 'noncyclic', families = Gumbel_Cop(marg1 = "LOGNO",
                                                                marg2 = "LOGLOG"))

  if(!any(selected(mod,merge = T) > p)) mod[5000]
  if(!any(selected(mod,merge = T) > p)) mod[10000]
  mstop_probes <- min(which(selected(mod,merge = T) > p   )) - 1 
  mod[mstop_probes]
  runtime = Sys.time() - time
  units(runtime) = "mins"
  
  mstop.mod<-  vector('list')
  mstop.mod$mstop <- mstop_probes
  mstop.mod$mu1 <- mod$mu1$mstop()
  mstop.mod$mu2 <- mod$mu2$mstop()
  mstop.mod$sigma1 <- mod$sigma1$mstop()
  mstop.mod$sigma2 <- mod$sigma2$mstop()
  mstop.mod$rho <- mod$rho$mstop()
  
  Coefficients <- coef(mod[mstop_probes])
  
  TP <- vector(length = 5, 'list')
  TP <- lapply(1:5,function(i){sum( gsub("^[^(]*\\(([^)]*)\\)","\\1", names(Coefficients[[i]])) %in% paste('x',which(TrueBeta[[i]]!=0),sep = ''))})
  names(TP) <- names(Coefficients)
  
  FP <- vector(length = 5, 'list')
  FP <- lapply(1:5,function(i){sum(!( gsub("^[^(]*\\(([^)]*)\\)","\\1", names(Coefficients[[i]])) %in% paste('x',which(TrueBeta[[i]]!=0),sep = '')))})
  names(FP) <- names(Coefficients)
  
  # Prediction error
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
  dLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
    
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
  pLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
    
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
  
  loss <- function (y, mu1, sigma1, mu2, sigma2, rho){
    -(log(dlnorm(x = y[, 1], meanlog = mu1, sdlog = sigma1)) +
        log(dLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2)) +
        log(rho + 1) - (1 + rho) * (log(plnorm(q = y[, 1], meanlog = mu1, sdlog = sigma1)) +
                                      log(pLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2))) -
        (2 + 1/rho) * log(plnorm(q = y[, 1], meanlog = mu1, sdlog = sigma1)^(-rho) +
                            pLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2)^(-rho) - 1))}
  LogLik_values <- loss(probes.test[,c(1,2)],mu1 = pred_mu1, mu2 = pred_mu2, 
                        sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho)
  
  LogLik <- sum(loss(probes.test[,c(1,2)],mu1 = pred_mu1, mu2 = pred_mu2, 
                     sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho), na.rm = T)
  
  es_copula <- rep(NA, 1000)
  
  for(pred_sim_run in 1:1000){
    
    pred_sample_copula <- matrix(NA, nrow = 2, ncol = 10000)
    
    Cop <- archmCopula(family = "gumbel", dim = 2, param = pred_rho[pred_sim_run])
    speclist1 <- list( meanlog = pred_mu1[pred_sim_run], sdlog = pred_sigma1[pred_sim_run] )
    speclist2 <- list( mu = pred_mu2[pred_sim_run], sigma = pred_sigma2[pred_sim_run], nu = 1, tau = 1)
    spec <- mvdc(copula = Cop, c("lnorm" , "GB2") , list(speclist1 , speclist2) )
    samp_cop <-  rMvdc(10000, spec)
    pred_sample_copula[1, ] <- samp_cop[ ,1]
    pred_sample_copula[2, ] <- samp_cop[ ,2]
    
    es_copula[pred_sim_run] <- es_sample(y = c(test$y1[pred_sim_run], test$y2[pred_sim_run]), dat = pred_sample_copula)
  }
  list(instance, Coefficients=Coefficients, TP=TP, FP=FP, mstop=mstop.mod, pred=pred, LogLik=LogLik, LogLik_values=LogLik_values, energy_score = es_copula,  runtime = runtime)
}



stability <- function(job, data, instance,q){

  train <- instance$train
  test <- instance$test
  weights_train <- instance$weights_train
  TrueBeta <- instance$TrueBeta
  seed = instance$seed
  train1 <- train[weights_train == 1,]
  
  set.seed(seed)
  form <- as.formula(( paste( "cbind(y1,y2) ~",  paste(c(paste("bbs(x", 1:p, ")", sep = "")), collapse = "+"))))
  
  time = Sys.time()
  mod = gamboostLSS(form, data = train1,
                    control = boost_control(mstop = 5*nrow(train), nu = 0.01),method = 'noncyclic', families = Gumbel_Cop(marg1 = "LOGNO",
                                                                                                                          marg2 = "LOGLOG"))
  
  if(p == 20) pfer = c(1,2.5,5)
  if(p == 1000) pfer = c(0.05, 0.1, 0.5)
  
  mstop.mod<-  vector('list')
  TP <- vector( 'list')
  FP <- vector( 'list')
  runtime_all <- vector('list')
  pred <- vector('list')
  LogLik <- vector('list')
  Coefficients_all <- vector('list')
  es_copula_ges <- vector("list")
  for(i in 1:length(pfer)){
    
    if(i == 1){
      s1 <- stabsel(mod, q =  q, PFER = pfer[1],assumption = "none")
      sel = names(s1$selected)  
      
    }else{
      cutoff <- min(1,((s1$q^2/(pfer[i]*s1$p)) + 1)/2)
      sel <- rownames(s1$phat)[s1$phat[,dim(s1$phat)[2]] > cutoff]
      time = Sys.time()
    }
    
    if(any(sel %in% (grep("(Intercept)",sel , value = T)))){
      sel <- sel[-which(sel %in% (grep("(Intercept)", sel, value = T)))]
    }
    
    
    if(any(grepl('.mu1', sel))){
      sel_mu1 <- gsub(".mu1","\\1",sel[grep('mu1',sel)])
      form_mu1 <- as.formula(paste("cbind(y1,y2)~",paste(sel_mu1,collapse = '+'), sep= ''))
    }else{
      sel_mu1 = 0
      form_mu1 <- as.formula(paste("cbind(y1,y2)~", 'bols(Intercept, intercept=F)'))
    }
    
    if(any(grepl('mu2', sel))){
      sel_mu2 <- gsub(".mu2","\\1",sel[grep('mu2', sel)])
      form_mu2 <- as.formula(paste("cbind(y1,y2)~",paste(sel_mu2,collapse = '+'),sep = ''))
    }else{
      sel_mu2 = 0
      form_mu2 <- as.formula(paste("cbind(y1,y2)~", 'bols(Intercept, intercept=F)'))
    }
    
    if(any(grepl('sigma1', sel))){
      sel_sigma1 <- gsub(".sigma1","\\1",sel[grep('sigma1', sel)])
      form_sigma1 <- as.formula(paste("cbind(y1,y2)~",paste(sel_sigma1,collapse = '+'), sep = ''))
      
    }else{
      sel_sigma1 = 0
      form_sigma1 <- as.formula(paste("cbind(y1,y2)~", 'bols(Intercept, intercept=F)'))
    }
    
    if(any(grepl('sigma2', sel))){
      sel_sigma2 <- gsub(".sigma2","\\1",sel[grep('sigma2', sel)])
      form_sigma2 <- as.formula(paste("cbind(y1,y2)~",paste(sel_sigma2,collapse = '+'),sep=''))
    }else{
      sel_sigma2 = 0
      form_sigma2 <- as.formula(paste("cbind(y1,y2)~", 'bols(Intercept, intercept=F)'))
    }
    
    if(any(grepl('rho', sel))){
      sel_rho <- gsub(".rho","\\1",sel[grep('rho', sel)])
      form_rho <- as.formula(paste("cbind(y1,y2)~",paste(sel_rho,collapse = '+'),sep = ''))
    }else{
      sel_rho = 0
      form_rho <- as.formula(paste("cbind(y1,y2)~", 'bols(Intercept, intercept=F)'))
    }
    # abfrage, ob 
    form <- list(mu1 = form_mu1,
                 mu2 = form_mu2,
                 sigma1 = form_sigma1,
                 sigma2 = form_sigma2,
                 rho = form_rho)
    
    mod = gamboostLSS(form, data = train,
                      control = boost_control(mstop = 5000, nu = 0.01,risk = 'oobag'), weights = weights_train, method = 'noncyclic', families = Gumbel_Cop(marg1 = "LOGNO",
                                                                                                                                                            marg2 = "LOGLOG"))
    
    MSTOP <- which.min(risk(mod,merge = T))
    oobag.risk <- risk(mod,merge = T)
    
    rm(mod) # removed the fist fitted model
    
    mod = gamboostLSS(form, data = train1, control = boost_control(mstop = MSTOP, nu = 0.01), method = 'noncyclic', families = Gumbel_Cop(marg1 = "LOGNO",marg2 = "LOGLOG"))
    
    runtime = Sys.time() - time
    units(runtime) = "mins"
    
    if(i == 1){
      runtime_all[[1]] <-  runtime 
    }else{
      runtime_all[[i]] <-  runtime
    }
    
    mstop.mod$mstop[[i]] <- MSTOP
    mstop.mod$mu1[[i]] <- mod$mu1$mstop()
    mstop.mod$mu2[[i]] <- mod$mu2$mstop()
    mstop.mod$sigma1[[i]] <- mod$sigma1$mstop()
    mstop.mod$sigma2[[i]] <- mod$sigma2$mstop()
    mstop.mod$rho[[i]] <- mod$rho$mstop()
    
    Coefficients_all[[i]] <- coef(mod)
    Coefficients <- coef(mod)
    
    TP[[i]] <- lapply(1:5,function(i){sum( gsub("\\s", "\\1", gsub("[bbs(),df=dfae]","\\1", names(Coefficients[[i]]))) %in% paste('x',which(TrueBeta[[i]]!=0),sep = ''))})
    names(TP[[i]]) <- names(Coefficients)
    
    FP[[i]] <- lapply(1:5,function(i){sum(gsub("\\s","\\1", gsub("[bbs(),df=dfae]","\\1", names(Coefficients[[i]]))) %in% paste('x',which(TrueBeta[[i]]==0),sep = ''))})
    names(FP[[i]]) <- names(Coefficients)
    
    #Prediction error
    pred_mu1    <- predict(mod$mu1, newdata = test, type = 'response')
    pred_sigma1 <- predict(mod$sigma1, newdata = test, type = 'response')
    pred_mu2    <- predict(mod$mu2, newdata = test, type = 'response')
    pred_sigma2 <- predict(mod$sigma2, newdata = test, type = 'response')
    pred_rho    <- predict(mod$rho, newdata = test, type = 'response')
    
    pred$mu1[[i]] <- pred_mu1
    pred$mu2[[i]] <- pred_mu2
    pred$sigma1[[i]] <- pred_sigma1
    pred$sigma2[[i]] <- pred_sigma2
    pred$rho[[i]] <- pred_rho
    
    ### helpers LOGLOG
    ### pdf 
    dLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
      
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
    pLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
      
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
    
    loss <- function (y, mu1, sigma1, mu2, sigma2, rho){
      -(log(dlnorm(x = y[, 1], meanlog = mu1, sdlog = sigma1)) +
          log(dLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2)) +
          log(rho + 1) - (1 + rho) * (log(plnorm(q = y[, 1], meanlog = mu1, sdlog = sigma1)) +
                                        log(pLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2))) -
          (2 + 1/rho) * log(plnorm(q = y[, 1], meanlog = mu1, sdlog = sigma1)^(-rho) +
                              pLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2)^(-rho) - 1))}
    
    LogLik[[i]] <- sum(loss(y = test[,c('y1','y2')], mu1 = pred_mu1, mu2 = pred_mu2,
                            sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho), na.rm = T)
    
    es_copula <- rep(NA, 1000)
    
    for(pred_sim_run in 1:1000){
      
      pred_sample_copula <- matrix(NA, nrow = 2, ncol = 10000)
      
      Cop <- archmCopula(family = "gumbel", dim = 2, param = pred_rho[pred_sim_run])
      speclist1 <- list( meanlog = pred_mu1[pred_sim_run], sdlog = pred_sigma1[pred_sim_run] )
      speclist2 <- list( mu = pred_mu2[pred_sim_run], sigma = pred_sigma2[pred_sim_run], nu = 1, tau = 1)
      spec <- mvdc(copula = Cop, c("lnorm" , "GB2") , list(speclist1 , speclist2) )
      samp_cop <-  rMvdc(10000, spec)
      pred_sample_copula[1, ] <- samp_cop[ ,1]
      pred_sample_copula[2, ] <- samp_cop[ ,2]
      
      es_copula[pred_sim_run] <- es_sample(y = c(test$y1[pred_sim_run], test$y2[pred_sim_run]), dat = pred_sample_copula)
      
    }
    es_copula_ges[[i]] <- es_copula
  }
  list(instance, Coefficients=Coefficients_all, TP=TP, FP=FP,  mstop=mstop.mod, pred=pred, LogLik=LogLik, energy_score = es_copula_ges, values = list(n=n, p=p, q=q, pfer=pfer), runtime = runtime_all)
}


deselection_change <- function(job, data, instance,tau_1){
  # data <- dat(seed, p = p,n = n)
  
  train <- instance$train
  test <- instance$test
  TrueBeta <- instance$TrueBeta
  weights_train <- instance$weights_train
  
  seed = instance$seed
  
  require(plyr)
  DeselectBoostLSS_5 <- function(object, data = NULL, fam, tau = NULL, method = c('attributable','cumulative')){
    data = attr(object,'data')
    
    mstop <- ifelse(inherits(object,'nc_mboostLSS'), mstop(object)[1],sum(mstop(object)))
    
    if (length(risk(object, merge = T)) > (mstop + 5)) stop("risk cannot contain more entries than mstop")
    
    parameter = names(object)
    which.response <- colnames(object[[1]]$response) == colnames(data)
    name.response <- colnames(data)[which.response]
    
    select <- selected(object, parameter = names(object))
    select1 <- selected(object, parameter = names(object))
    if(inherits(object,'glmboostLSS')){
      select[[1]] <- select[[1]]-1  
      select[[2]] <- select[[2]]-1 
      select[[3]] <- select[[3]]-1
      select[[4]] <- select[[4]]-1
      select[[5]] <- select[[5]]-1
      select[[1]] <- select[[1]][select[[1]] !=0]
      select[[2]] <- select[[2]][select[[2]] !=0]
      select[[3]] <- select[[3]][select[[3]] !=0]
      select[[4]] <- select[[4]][select[[4]] !=0]
      select[[5]] <- select[[5]][select[[5]] !=0]
      
      RiskRed <- risk(object,merge = T)
      totalRiskRed <- (risk(object, merge=T)[1]-risk(object, merge=T)[length(risk(object, merge=T))])
      diffRiskRed_all <-sapply(5:(length(RiskRed)-1),function(k){RiskRed[k]-RiskRed[k+1]})
      diffRiskRed <- vector("list")
      diffRiskRed[[1]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[1]]][select1[[1]]-1 != 0]
      diffRiskRed[[2]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[2]]][select1[[2]]-1 != 0]
      diffRiskRed[[3]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[3]]][select1[[3]]-1 != 0]
      diffRiskRed[[4]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[4]]][select1[[4]]-1 != 0]
      diffRiskRed[[5]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[5]]][select1[[5]]-1 != 0]
      
    }else{
      RiskRed <- risk(object,merge = T)
      diffRiskRed_all <-sapply(1:(length(RiskRed)-1),function(k){RiskRed[k]-RiskRed[k+1]})
      diffRiskRed_all <- diffRiskRed_all[-c(1:4)]
      names(diffRiskRed_all) <- names(risk(object,merge= T)[-c(1:5)])
      
      diffRiskRed <- vector("list")
      diffRiskRed[[1]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[1]]]
      diffRiskRed[[2]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[2]]]
      diffRiskRed[[3]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[3]]]
      diffRiskRed[[4]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[4]]]
      diffRiskRed[[5]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[5]]]
    }
    
    if(is.null(coef(object)[[1]])) select[[1]]  = 0 else  select[[1]] = select[[1]] 
    if(is.null(coef(object)[[2]])) select[[2]]  = 0 else  select[[2]] = select[[2]] 
    if(is.null(coef(object)[[3]])) select[[3]]  = 0 else  select[[3]] = select[[3]] 
    if(is.null(coef(object)[[4]])) select[[4]]  = 0 else  select[[4]] = select[[4]] 
    if(is.null(coef(object)[[5]])) select[[5]]  = 0 else  select[[5]] = select[[5]] 
    
    Var = lapply(1:length(parameter), function(i){count(select[[i]])[2]})
    Risk.Var <- lapply(1:5,function(i){sapply(1:dim(Var[[i]])[1],function(j){sum(diffRiskRed[[i]][which(count(select[[i]])[[1]][j] == select[[i]])])})})
    
    w.parameter <- c(rep(parameter[1],length(count(select[[1]])[[1]])),rep(parameter[2],length(count(select[[2]])[[1]])),rep(parameter[3],length(count(select[[3]])[[1]])),rep(parameter[4],length(count(select[[4]])[[1]])),rep(parameter[5],length(count(select[[5]])[[1]])))
    n.parameter <- c(names(object[[1]]$coef()),names(object[[2]]$coef()),names(object[[3]]$coef()),names(object[[4]]$coef()),names(object[[5]]$coef()))
    if('(Intercept)' %in% n.parameter)  n.parameter<-n.parameter[-which(n.parameter == '(Intercept)')]
    
    Risk.order <- data.frame(w.parameter[!is.na(unlist(Risk.Var))], n.parameter, unlist(Risk.Var)[!is.na(unlist(Risk.Var))])
    colnames(Risk.order) <- c('parameter', 'VarName', 'Risk')
    Risk.order <- Risk.order[order(Risk.order$Risk),]
    Risk.order$CumRisk <- cumsum(Risk.order$Risk)
    
    perc <- ifelse(is.null(tau), 0.01, tau) 
    percRiskRed <-  (risk(object, merge=T)[1]-risk(object, merge=T)[length(risk(object, merge=T))]) * perc
    if(method[1] == 'attributable'){RiskRedOver <- Risk.order[which(Risk.order$Risk > percRiskRed),]
    }else{RiskRedOver <- Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
    
    if(is.na(parameter[2])) parameter[2]<-0
    if(is.na(parameter[1])) parameter[1]<-0
    if(is.na(parameter[3])) parameter[3]<-0
    if(is.na(parameter[4])) parameter[4]<-0
    if(is.na(parameter[5])) parameter[5]<-0
    
    if(length(name.response) > 1) name.response <- paste("cbind(",paste(name.response, collapse = ","),")")
    
    if(any(RiskRedOver$parameter == parameter[1])){
      form1 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[1]], collapse= "+")))
    }else{  form1 <-  as.formula(paste(name.response, '~bols(Intercept, intercept=F)'))}
    
    if(any(RiskRedOver$parameter == parameter[2])){
      form2 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[2]], collapse= "+"))) 
    }else{  form2 <-  as.formula(paste(name.response, '~bols(Intercept, intercept=F)'))}
    
    if(any(RiskRedOver$parameter == parameter[3])){
      form3 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[3]], collapse= "+"))) 
    }else{  form3 <-  as.formula(paste(name.response, '~bols(Intercept, intercept=F)'))}
    
    if(any(RiskRedOver$parameter == parameter[4])){
      form4 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[4]], collapse= "+"))) 
    }else{  form4 <-  as.formula(paste(name.response, '~bols(Intercept, intercept=F)'))}
    
    if(any(RiskRedOver$parameter == parameter[5])){
      form5 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[5]], collapse= "+"))) 
    }else{  form5 <-  as.formula(paste(name.response, '~bols(Intercept, intercept=F)'))}
    
    formula <- list(form1,form2,form3,form4,form5) 
    names(formula)<- names(object)
    dfbase <- environment(environment(environment(object[[1]][["fitted"]])[["RET"]][["baselearner"]][[1]][["model.frame"]])[["ret"]][["model.frame"]])[["df"]]
    
    
    if(inherits(object,"nc_mboostLSS")){
      if(inherits(object,"glmboostLSS")){
        model_after = glmboostLSS(formula, data = data, families = fam, method = 'noncyclic', weights = model.weights(object), control = boost_control(mstop = mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu),as.numeric(object[[4]]$control$nu),as.numeric(object[[5]]$control$nu))))
      }else{
        model_after = gamboostLSS(formula, data = data, families = fam, method = 'noncyclic', weights = model.weights(object), control = boost_control(mstop =  mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu),as.numeric(object[[4]]$control$nu),as.numeric(object[[5]]$control$nu))))
      }
    }else{
      if(inherits(object,"glmboostLSS")){
        model_after = glmboostLSS(formula, data = data, families = fam, method = 'cyclic', weights = model.weights(object), control = boost_control(mstop = mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu),as.numeric(object[[4]]$control$nu),as.numeric(object[[5]]$control$nu))))
      }else{
        model_after = gamboostLSS(formula, data = data, families = fam, method = 'cyclic', weights = model.weights(object), control = boost_control(mstop =  mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu),as.numeric(object[[4]]$control$nu),as.numeric(object[[5]]$control$nu))))
      }
    }
    
    out <- model_after
    Coef <- coef(model_after)
    deselect_para <-  list(coef = Coef, tau = perc, deselectmethod = method[1])
    out <- append(x = out, values = deselect_para)
    
    class(out) <- c(class(out))
    
    return(out)
  }
  
  form <- as.formula(( paste( "cbind(y1,y2) ~",  paste(c(paste("bbs(x", 1:p, ")", sep = "")), collapse = "+"))))
  
  set.seed(seed)
  
  time = Sys.time()
  mod = gamboostLSS(form, data = train, control = boost_control(mstop = 5000, nu = 0.01, risk = "oobag"), weights = weights_train, method = 'noncyclic', families = Gumbel_Cop(marg1 = "LOGNO",
                                                                                                                                                                               marg2 = "LOGLOG"))
  MSTOP <- which.min(risk(mod,merge = T))
  
  oobag.risk <- risk(mod,merge = T)
  
  rm(mod) # removed the fist fitted model
  train <- train[weights_train == 1, ]
  
  mod = gamboostLSS(form, data = train, control = boost_control(mstop = MSTOP, nu = 0.01), method = 'noncyclic', families = Gumbel_Cop(marg1 = "LOGNO",marg2 = "LOGLOG"))
  
  time = Sys.time()
  mod_desel <- DeselectBoostLSS_5(mod, tau = tau_1, fam = Gumbel_Cop(marg1 = "LOGNO",
                                                                     marg2 = "LOGLOG"))
  
  runtime = Sys.time() - time
  units(runtime) = "mins"
  
  mstop.mod<-  vector('list')
  mstop.mod$mstop <- MSTOP
  mstop.mod$mu1 <- mod$mu1$mstop()
  mstop.mod$mu2 <- mod$mu2$mstop()
  mstop.mod$sigma1 <- mod$sigma1$mstop()
  mstop.mod$sigma2 <- mod$sigma2$mstop()
  mstop.mod$rho <- mod$rho$mstop()
  
  Coefficients <- mod_desel$coef
  
  TP <- vector(length = 5, 'list')
  TP <- lapply(1:5,function(i){sum( gsub("^[^(]*\\(([^)]*)\\)","\\1", names(Coefficients[[i]])) %in% paste('x',which(TrueBeta[[i]]!=0),sep = ''))})
  names(TP) <- names(Coefficients)
  
  FP <- vector(length = 5, 'list')
  FP <- lapply(1:5,function(i){sum(!( gsub("^[^(]*\\(([^)]*)\\)","\\1", names(Coefficients[[i]])) %in% paste('x',which(TrueBeta[[i]]!=0),sep = '')))})
  names(FP) <- names(Coefficients)
  
  #Prediction error
  pred_mu1    <- predict(mod_desel$mu1, newdata = test, type = 'response')
  pred_sigma1 <- predict(mod_desel$sigma1, newdata = test, type = 'response')
  pred_mu2    <- predict(mod_desel$mu2, newdata = test, type = 'response')
  pred_sigma2 <- predict(mod_desel$sigma2, newdata = test, type = 'response')
  pred_rho    <- predict(mod_desel$rho, newdata = test, type = 'response')
  
  pred <- vector('list')
  pred$mu1 <- pred_mu1
  pred$mu2 <- pred_mu2
  pred$sigma1 <- pred_sigma1
  pred$sigma2 <- pred_sigma2
  pred$rho <- pred_rho
  
  ### helpers LOGLOG
  ### pdf 
  dLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
    
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
  pLogLogistic <- function (x, mu = 1, sigma = 1, log = FALSE){
    
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
  
  loss <- function (y, mu1, sigma1, mu2, sigma2, rho){
    -(log(dlnorm(x = y[, 1], meanlog = mu1, sdlog = sigma1)) +
        log(dLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2)) +
        log(rho + 1) - (1 + rho) * (log(plnorm(q = y[, 1], meanlog = mu1, sdlog = sigma1)) +
                                      log(pLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2))) -
        (2 + 1/rho) * log(plnorm(q = y[, 1], meanlog = mu1, sdlog = sigma1)^(-rho) +
                            pLogLogistic(x = y[, 2], mu = mu2, sigma = sigma2)^(-rho) - 1))}
  LogLik <- sum(loss(y = test[,c('y1','y2')], mu1 = pred_mu1, mu2 = pred_mu2,
                     sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho), na.rm = T)
  
  es_copula <- rep(NA, 1000)
  
  for(pred_sim_run in 1:1000){
    
    pred_sample_copula <- matrix(NA, nrow = 2, ncol = 10000)
    
    Cop <- archmCopula(family = "gumbel", dim = 2, param = pred_rho[pred_sim_run])
    speclist1 <- list( meanlog = pred_mu1[pred_sim_run], sdlog = pred_sigma1[pred_sim_run] )
    speclist2 <- list( mu = pred_mu2[pred_sim_run], sigma = pred_sigma2[pred_sim_run], nu = 1, tau = 1)
    spec <- mvdc(copula = Cop, c("lnorm" , "GB2") , list(speclist1 , speclist2) )
    samp_cop <-  rMvdc(10000, spec)
    pred_sample_copula[1, ] <- samp_cop[ ,1]
    pred_sample_copula[2, ] <- samp_cop[ ,2]
    
    es_copula[pred_sim_run] <- es_sample(y = c(test$y1[pred_sim_run], test$y2[pred_sim_run]), dat = pred_sample_copula)
    
  } 
  list(instance, Coefficients=Coefficients, TP=TP, FP=FP, mstop=mstop.mod, pred=pred, LogLik=LogLik, energy_score = es_copula, runtime = runtime)
}


addProblem(reg = reg, "data123", seed = 124,data = list(p=p,n=n), fun = dat)
addAlgorithm(reg = reg,'classic', fun = classic)
addAlgorithm(reg = reg, 'probing', fun = probing)
addAlgorithm(reg = reg, 'deselection_change',fun = deselection_change)
addAlgorithm(reg = reg,'stability', fun = stability)

algo.design <- list()
algo.design$classic <- data.frame()
algo.design$probing <- data.frame()
algo.design$deselection_change <- expand.grid(tau_1 = TAU_1)
algo.design$stability <- expand.grid(q=Q)


addExperiments(reg = reg, algo.designs = algo.design, 
               repls = REPLICATIONS)

summarizeExperiments()
submitJobs()
getStatus()

