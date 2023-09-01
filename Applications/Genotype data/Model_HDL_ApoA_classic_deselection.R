################################################################
################################################################
# ---------- Application genotype data HDL and ApoA  --------- #

library("Rcpp")
library("gamboostLSS")
library("gamlss.dist")
library("copula")

source("defs.R")
source("Copulas/Copula_Gaussian.R")
source("Copulas/Copula_Clayton.R")
source("Copulas/Copula_Gumbel.R")

sim <- function(seed){

tau_1 = 0.01
load("data_HDL.RData")
data_HDL <- data[[1]]$HDL

load("data_ApoA.RData")
dataset <- data[[1]]
weights <- data[[2]]
dataset <- cbind(HDL = data_HDL,dataset)

load("data_HDL_ApoA_test.RData")

train = dataset
test = test 
weights_train=weights

set.seed(seed)

time_classic = Sys.time()
mod = glmboostLSS(cbind(HDL,ApoA)~., data = train, control = boost_control(mstop = 2000, nu = 0.01, risk = "oobag", trace = T), weights = weights_train, method = 'noncyclic', families =  Gauss_Cop(marg1 = "LOGNO", marg2 = "LOGNO", stabilization = "MAD"))
MSTOP <- which.min(risk(mod,merge = T))# seed = 1, MSTOP = 537

if(MSTOP >= 990){
  mod[2000]
}
MSTOP <- which.min(risk(mod,merge = T))

if(MSTOP >= 1990){
  mod[3000]
}
MSTOP <- which.min(risk(mod,merge = T))

if(MSTOP >= 2990){
  mod[4000]
}
MSTOP <- which.min(risk(mod,merge = T))

if(MSTOP >= 3990){
  mod[5000]
}
MSTOP <- which.min(risk(mod,merge = T))

if(MSTOP >= 4990){
  mod[6000]
}
MSTOP <- which.min(risk(mod,merge = T))


if(MSTOP >= 5990){
  mod[7000]
}
MSTOP <- which.min(risk(mod,merge = T))

if(MSTOP >= 6990){
  mod[8000]
}
MSTOP <- which.min(risk(mod,merge = T))


if(MSTOP >= 7990){
  mod[9000]
}
MSTOP <- which.min(risk(mod,merge = T))


if(MSTOP >= 8990){
  mod[10000]
}
MSTOP <- which.min(risk(mod,merge = T))

if(MSTOP >= 9990){
  mod[12000]
}
MSTOP <- which.min(risk(mod,merge = T))

if(MSTOP >= 11990){
  mod[14000]
}
MSTOP <- which.min(risk(mod,merge = T))

oobag.risk <- risk(mod,merge = T)

rm(mod) # removed the fist fitted model
train <- train[weights_train == 1, ]

mod = glmboostLSS(cbind(HDL,ApoA)~., data = train, control = boost_control(mstop = MSTOP, nu = 0.01, trace = T), method = 'noncyclic', families = Gauss_Cop(marg1 = "LOGNO", marg2 = "LOGNO", stabilization = "MAD"))

runtime_classic = Sys.time() - time_classic
units(runtime_classic) = "mins"

mstop.mod<-  vector('list')
mstop.mod$mstop <- MSTOP
mstop.mod$mu1 <- mod$mu1$mstop()
mstop.mod$mu2 <- mod$mu2$mstop()
mstop.mod$sigma1 <- mod$sigma1$mstop()
mstop.mod$sigma2 <- mod$sigma2$mstop()
mstop.mod$rho <- mod$rho$mstop()

Coefficients <- coef(mod)

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
  -(log(dlnorm(x = y[, 1], meanlog = mu1, sdlog = sigma1)) + log(dlnorm(x = y[,  
                                                                            2], meanlog = mu2, sdlog = sigma2)) - 0.5 * log(1 - rho^2) +  
      (rho/(1 - rho^2)) * qnorm(plnorm(q = y[, 1], meanlog = mu1,  
                                       sdlog = sigma1)) * qnorm(plnorm(q = y[, 2], meanlog = mu2,  
                                                                       sdlog = sigma2)) - (rho^2/(2 * (1 - rho^2))) * (qnorm(plnorm(q = y[, 1], meanlog = mu1, sdlog = sigma1))^2 + qnorm(plnorm(q = y[,   2], meanlog = mu2, sdlog = sigma2))^2)) 
  }


LogLik_values <- loss(y = test[,c('HDL','ApoA')], mu1 = pred_mu1, mu2 = pred_mu2,
                      sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho)
LogLik <- sum(loss(y = test[,c('HDL','ApoA')], mu1 = pred_mu1, mu2 = pred_mu2,
                   sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho),na.rm=T)


res <- list(mod = mod, Coefficients=Coefficients, mstop=mstop.mod, pred=pred, LogLik_values=LogLik_values, 
            runtime_classic = runtime_classic)

save(res, file = paste0("sims_HDL_ApoA__classic_1_MAD_",seed,".RData"))
rm(res, train, weights_train)
gc(reset=T)


###############################
###############################
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
  
  RiskRed <- risk(object,merge = T)
  totalRiskRed <- (risk(object, merge=T)[1]-risk(object, merge=T)[length(risk(object, merge=T))])
  diffRiskRed_all <-sapply(1:(length(RiskRed)-1),function(k){RiskRed[k]-RiskRed[k+1]})
  diffRiskRed <- vector("list")
  diffRiskRed[[1]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[1]]]
  diffRiskRed[[2]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[2]]]
  diffRiskRed[[3]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[3]]]
  diffRiskRed[[4]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[4]]]
  diffRiskRed[[5]] <- diffRiskRed_all[names(diffRiskRed_all)==parameter[[5]]]
  
  
  Var = lapply(1:length(parameter), function(i){count(select[[i]])[2]})

  Risk.Var <- lapply(1:length(parameter),function(i){sapply(1:dim(Var[[i]])[1],function(j){sum(diffRiskRed[[i]][which(count(select[[i]])[[1]][j] == select[[i]])])})})

  w.parameter <- c(rep(parameter[1],length(count(select[[1]])[[1]])),rep(parameter[2],length(count(select[[2]])[[1]])),rep(parameter[3],length(count(select[[3]])[[1]])),rep(parameter[4],length(count(select[[4]])[[1]])),rep(parameter[5],length(count(select[[5]])[[1]])))
  n.parameter <- c(names(object[[1]]$coef()),names(object[[2]]$coef()),names(object[[3]]$coef()),names(object[[4]]$coef()),names(object[[5]]$coef()))
  if('(Intercept)' %in% n.parameter) {
    w.parameter <- w.parameter[-which(n.parameter == '(Intercept)')]
    Risk.Var <- unlist(Risk.Var)[-which(n.parameter == '(Intercept)')]
    n.parameter<-n.parameter[-which(n.parameter == '(Intercept)')]
  }
  
  Risk.order <- data.frame(w.parameter, n.parameter, unlist(Risk.Var))
  colnames(Risk.order) <- c('parameter', 'VarName', 'Risk')
  Risk.order <- Risk.order[order(Risk.order$Risk),]
  Risk.order$CumRisk <- cumsum(Risk.order$Risk)
  
  perc <- ifelse(is.null(tau), 0.01, tau)
  percRiskRed <-  sum(Risk.order$Risk) * perc
  if(method[1] == 'attributable'){RiskRedOver <- Risk.order[which(Risk.order$Risk > percRiskRed),]
  }else{RiskRedOver <- Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
  
  if(is.na(parameter[2])) parameter[2]<-0
  if(is.na(parameter[1])) parameter[1]<-0
  if(is.na(parameter[3])) parameter[3]<-0
  if(is.na(parameter[4])) parameter[4]<-0
  if(is.na(parameter[5])) parameter[5]<-0
  
  if(length(name.response) > 1) name.response <- paste("cbind(",paste(name.response, collapse = ","),")")
  
  if(any(RiskRedOver$parameter == parameter[1])){
    # if(length(name.response) >1) name.response <- paste("cbind(",paste(name.response, collapse = ","),")")
    form1 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[1]], collapse= "+")))
  }else{  form1 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[2])){
    form2 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[2]], collapse= "+")))
  }else{  form2 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[3])){
    form3 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[3]], collapse= "+")))
  }else{  form3 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[4])){
    form4 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[4]], collapse= "+")))
  }else{  form4 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[5])){
    form5 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[5]], collapse= "+")))
  }else{  form5 <-  as.formula(paste(name.response, '~1'))}
  
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



time_desel = Sys.time()

mod_desel <- DeselectBoostLSS_5(mod, tau = tau_1, fam = Gauss_Cop(marg1 = "LOGNO", marg2 = "LOGNO", stabilization = "MAD"))

runtime_desel = Sys.time() - time_desel
units(runtime_desel) = "mins"


Coefficients_desel <- mod_desel$coef

#Prediction error
pred_mu1_desel    <- predict(mod_desel$mu1, newdata = test, type = 'response')
pred_sigma1_desel <- predict(mod_desel$sigma1, newdata = test, type = 'response')
pred_mu2_desel    <- predict(mod_desel$mu2, newdata = test, type = 'response')
pred_sigma2_desel <- predict(mod_desel$sigma2, newdata = test, type = 'response')
pred_rho_desel    <- predict(mod_desel$rho, newdata = test, type = 'response')

pred_desel <- vector('list')
pred_desel$mu1 <- pred_mu1_desel
pred_desel$mu2 <- pred_mu2_desel
pred_desel$sigma1 <- pred_sigma1_desel
pred_desel$sigma2 <- pred_sigma2_desel
pred_desel$rho <- pred_rho_desel


LogLik_desel <- sum(loss(y = test[,c('HDL','ApoA')], mu1 = pred_mu1_desel, mu2 = pred_mu2_desel,
                         sigma1 = pred_sigma1_desel, sigma2 = pred_sigma2_desel, rho = pred_rho_desel), na.rm = T)



res <- list(Coefficients=Coefficients_desel, mstop=mstop.mod, pred=pred_desel, LogLik=LogLik_desel,
            runtime= runtime_desel)
save(res, file = paste0("sims_hdl_apoa_desel",seed,".RData"))



return(list(Coefficients=Coefficients, Coefficients_desel = Coefficients_desel, mstop=mstop.mod, pred=pred, pred_desel=pred_desel, pred_plus=pred_plus, LogLik_values=LogLik_values, LogLik=LogLik, LogLik_desel = LogLik_desel, runtime_classic = runtime_classic, runtime_desel = runtime_desel))
}


results = mclapply(1, sim, mc.cores = 1, mc.preschedule = FALSE)

save(results, file="HDL_ApoA_models.RData")

