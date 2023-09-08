################################################################
################################################################
# ------------ Analysis of fetal ultrasound data ------------- #

library("BatchJobs")
library(Rcpp)
library(gamboostLSS)
library(gamlss.dist)
library(copula)


source("defs.R")
source("Copulas/Copula_Gaussian.R")
source("Copulas/Copula_Clayton.R")
source("Copulas/Copula_Gumbel.R")




### load data

load("data_nmv.RData")


#View(data_nmv)
any(is.na(data_nmv))
dim(data_nmv)
names(data_nmv)
summary(data_nmv)

sapply(data_nmv, class)


# building a subset only including US measurements with a scan-to-delivery interval of 7 days
# no Intrauterine fetal deaths (stillbirths) and fetuses with structural anomalies.

data_subset <- data_nmv %>% 
  filter(INTERVAL < 8) %>% 
  filter(FEHLBILD_new == "nein") %>%
  filter(TOTGEBURT_new == "nein")

dim(data_subset)  

#View(data_subset)

# 6103 observations remaining

summary(data_subset)

names(data_subset)

data_subset <- data_subset[, c("GEBLAENGEK", "GEWICHT", 
                               "ABDOMEN.SAGITTALDURCHMESSER", "ABDOMEN.TRANSVERSALDURCHMESSE", "ABDOMENUMFANG",
                               "BIPARIETALER.DURCHMESSER", "FRONTOOKZIPITALER.DURCHMESSER", "KOPFUMFANG",
                               "FEMURLAENGE",
                               "ENDGEWICHT", "MUTTERGROESSE", "bmi_mutter","GDiabetis", "GRAVIDA", "PARA",
                               "GestAge_Days_US", "GESCHLECHT")]

data_subset
dim(data_subset)
names(data_subset)

### scaled Birth weight
data_subset["GEWICHT_sca"] <- data_subset["GEWICHT"] / 1000
summary(data_subset)

### interaction terms of the sonograpic variables

data_subset["ABSAG_ABTRAN"] <- data_subset["ABDOMEN.SAGITTALDURCHMESSER"] * data_subset["ABDOMEN.TRANSVERSALDURCHMESSE"]
data_subset["ABSAG_ABUMF"] <- data_subset["ABDOMEN.SAGITTALDURCHMESSER"] * data_subset["ABDOMENUMFANG"]
data_subset["ABSAG_BIPADU"] <- data_subset["ABDOMEN.SAGITTALDURCHMESSER"] * data_subset["BIPARIETALER.DURCHMESSER"]
data_subset["ABSAG_FRODU"] <- data_subset["ABDOMEN.SAGITTALDURCHMESSER"] * data_subset["FRONTOOKZIPITALER.DURCHMESSER"]
data_subset["ABSAG_KOPUM"] <- data_subset["ABDOMEN.SAGITTALDURCHMESSER"] * data_subset["KOPFUMFANG"]
data_subset["ABSAG_FEMLA"] <- data_subset["ABDOMEN.SAGITTALDURCHMESSER"] * data_subset["FEMURLAENGE"]

data_subset["ABTRAN_ABDUM"] <- data_subset["ABDOMEN.TRANSVERSALDURCHMESSE"] * data_subset["ABDOMENUMFANG"]
data_subset["ABTRAN_BIPADU"] <- data_subset["ABDOMEN.TRANSVERSALDURCHMESSE"] * data_subset["BIPARIETALER.DURCHMESSER"]
data_subset["ABTRAN_FRONDU"] <- data_subset["ABDOMEN.TRANSVERSALDURCHMESSE"] * data_subset["FRONTOOKZIPITALER.DURCHMESSER"]
data_subset["ABTRAN_KOPUM"] <- data_subset["ABDOMEN.TRANSVERSALDURCHMESSE"] * data_subset["KOPFUMFANG"]
data_subset["ABTRAN_FEMLA"] <- data_subset["ABDOMEN.TRANSVERSALDURCHMESSE"] * data_subset["FEMURLAENGE"]

data_subset["ABDUM_BIPADU"] <- data_subset["ABDOMENUMFANG"] * data_subset["BIPARIETALER.DURCHMESSER"]
data_subset["ABDUM_FRONDU"] <- data_subset["ABDOMENUMFANG"] * data_subset["FRONTOOKZIPITALER.DURCHMESSER"]
data_subset["ABD_KOPUM"] <- data_subset["ABDOMENUMFANG"] * data_subset["KOPFUMFANG"]
data_subset["ABDUM_FEMLA"] <- data_subset["ABDOMENUMFANG"] * data_subset["FEMURLAENGE"]

data_subset["BIPADU_FRODU"] <- data_subset["BIPARIETALER.DURCHMESSER"] * data_subset["FRONTOOKZIPITALER.DURCHMESSER"]
data_subset["BIPADU_KOPUM"] <- data_subset["BIPARIETALER.DURCHMESSER"] * data_subset["KOPFUMFANG"]
data_subset["BIPADU_FEMLA"] <- data_subset["BIPARIETALER.DURCHMESSER"] * data_subset["FEMURLAENGE"]

data_subset["FRONDU_KOPUM"] <- data_subset["FRONTOOKZIPITALER.DURCHMESSER"] * data_subset["KOPFUMFANG"]
data_subset["FRONDU_FEMLA"] <- data_subset["FRONTOOKZIPITALER.DURCHMESSER"] * data_subset["FEMURLAENGE"]

data_subset["KOPUM_FEMLA"] <- data_subset["KOPFUMFANG"] * data_subset["FEMURLAENGE"]


dim(data_subset)
names(data_subset)
summary(data_subset)
sapply(data_subset, class)


### histograms of the responses
hist(data_subset$GEWICHT, breaks = 50)
hist(data_subset$GEWICHT_sca, breaks = 50)

hist(data_subset$GEBLAENGEK, breaks = 50)





################################################################################
########################## Boosting and CV #####################################
################################################################################

# split the data set into training and test 

set.seed(1)

weights <- sample(c(rep(1, 4103), rep(0, 2000)))
data_subset_train <- data_subset[weights == 1, ]

### formula
names(data_subset)

mu1_eq <- cbind(GEBLAENGEK, GEWICHT_sca) ~ bols(GDiabetis) + 
  bols(GESCHLECHT) + 
  bbs(ABDOMEN.SAGITTALDURCHMESSER) + 
  bbs(ABDOMEN.TRANSVERSALDURCHMESSE) + 
  bbs(ABDOMENUMFANG) + 
  bbs(BIPARIETALER.DURCHMESSER) + 
  bbs(FRONTOOKZIPITALER.DURCHMESSER) + 
  bbs(KOPFUMFANG) + 
  bbs(FEMURLAENGE) +
  bbs(ABSAG_ABTRAN) +
  bbs(ABSAG_ABUMF) +
  bbs(ABSAG_BIPADU) +
  bbs(ABSAG_FRODU) +
  bbs(ABSAG_KOPUM) +
  bbs(ABSAG_FEMLA) +
  bbs(ABTRAN_ABDUM) +
  bbs(ABTRAN_BIPADU) + 
  bbs(ABTRAN_FRONDU) +
  bbs(ABTRAN_KOPUM) +
  bbs(ABTRAN_FEMLA) + 
  bbs(ABDUM_BIPADU) +
  bbs(ABDUM_FRONDU) +
  bbs(ABD_KOPUM) +
  bbs(ABDUM_FEMLA) +
  bbs(BIPADU_FRODU) +
  bbs(BIPADU_KOPUM) +
  bbs(BIPADU_FEMLA) +
  bbs(FRONDU_KOPUM) +
  bbs(FRONDU_FEMLA) +
  bbs(KOPUM_FEMLA) +
  bbs(ENDGEWICHT) + 
  bbs(MUTTERGROESSE) + 
  bbs(bmi_mutter) + 
  bbs(GRAVIDA) +
  bbs(PARA) +
  bbs(GestAge_Days_US) 

sigma1_eq <- mu2_eq <- sigma2_eq <- rho_eq <- mu1_eq

eq_list <- list(mu1 = mu1_eq,
                sigma1 = sigma1_eq,
                mu2 = mu2_eq,
                sigma2 = sigma2_eq,
                rho = rho_eq)



####################### Gauss: tanh link function ##############################
gauss_tanh_boost_fit <- gamboostLSS(eq_list,
                                    data = data_subset,
                                    weights = weights,
                                    families = Gauss_Cop(marg1 = "LOGLOG",
                                                         marg2 = "LOGLOG",
                                                         sigma1 = sd(data_subset_train$GEBLAENGEK),
                                                         sigma2 = sd(data_subset_train$GEWICHT_sca),
                                                         stabilization = "L2"),
                                    method = "noncyclic",
                                    control = boost_control(mstop = 8000,
                                                            nu = 0.01,
                                                            trace = T))

coef(gauss_tanh_boost_fit)
names(coef(gauss_tanh_boost_fit)$mu1)
names(coef(gauss_tanh_boost_fit)$sigma1)
names(coef(gauss_tanh_boost_fit)$mu2)
names(coef(gauss_tanh_boost_fit)$sigma2)
names(coef(gauss_tanh_boost_fit)$rho)



# tuning via noncyclic CV
set.seed(2)

folds <- cv(model.weights(gauss_tanh_boost_fit), type = "kfold", B = 10)
head(folds)
colSums(folds)

CV_mod <- cvrisk(gauss_tanh_boost_fit, folds = folds)

mstop(CV_mod)
plot(CV_mod)

# run final model on optimal mstop
gauss_tanh_boost_fit_final <- gamboostLSS(eq_list,
                                          data = data_subset,
                                          families = Gauss_Cop(marg1 = "LOGLOG",
                                                               marg2 = "LOGLOG",
                                                               sigma1 = sd(data_subset_train$GEBLAENGEK),
                                                               sigma2 = sd(data_subset_train$GEWICHT_sca),
                                                               stabilization = "L2"),
                                          weights = weights,
                                          method = "noncyclic",
                                          control = boost_control(mstop = mstop(CV_mod),
                                                                  nu = 0.01,
                                                                  trace = T))

object = gauss_tanh_boost_fit_final


pred_mu1    <- predict(object$mu1, newdata = test, type = 'response')
pred_sigma1 <- predict(object$sigma1, newdata = test, type = 'response')
pred_mu2    <- predict(object$mu2, newdata = test, type = 'response')
pred_sigma2 <- predict(object$sigma2, newdata = test, type = 'response')
pred_rho    <- predict(object$rho, newdata = test, type = 'response')

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

loss <-  function(y, mu1, sigma1, mu2, sigma2, rho){-(log(dLogLogistic(x = y[, 1], mu = mu1, sigma = sigma1)) + log(dLogLogistic(x = y[,2], mu = mu2, sigma = sigma2)) - 0.5 * log(1 - tanh(rho)^2) +  
                                                        (tanh(rho)/(1 - tanh(rho)^2)) * qnorm(pLogLogistic(x = y[, 1],mu = mu1, sigma = sigma1)) * qnorm(pLogLogistic(x = y[,2], mu = mu2, sigma = sigma2)) - (tanh(rho)^2/(2 * (1 - tanh(rho)^2))) * (qnorm(pLogLogistic(x = y[, 1], mu = mu1,  
                                                                                                                                                                                                                                                                          sigma = sigma1))^2 + qnorm(pLogLogistic(x = y[, 2], mu = mu2,  
                                                                                                                                                                                                                                                                                                                  sigma = sigma2))^2)) }
LogLik <- sum(loss(y = test[,c("GEBLAENGEK", "GEWICHT_sca")], mu1 = pred_mu1, mu2 = pred_mu2,
                         sigma1 = pred_sigma1, sigma2 = pred_sigma2, rho = pred_rho), na.rm = T)


########################################
############ Deselection ###############

library(plyr)
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
  percRiskRed <-  sum(Risk.order$Risk) * perc
  if(method[1] == 'attributable'){
    RiskRedOver <- Risk.order[which(Risk.order$Risk > percRiskRed),]
  }else{RiskRedOver <- Risk.order[which(Risk.order$CumRisk > percRiskRed),]}
  
  if(is.na(parameter[2])) parameter[2]<-0
  if(is.na(parameter[1])) parameter[1]<-0
  if(is.na(parameter[3])) parameter[3]<-0
  if(is.na(parameter[4])) parameter[4]<-0
  if(is.na(parameter[5])) parameter[5]<-0
  
  if(length(name.response) > 1) name.response <- paste("cbind(",paste(name.response, collapse = ","),")")
  
  data$Intercept <- 1
  if(any(RiskRedOver$parameter == parameter[1])){
    # if(length(name.response) >1) name.response <- paste("cbind(",paste(name.response, collapse = ","),")")
    form1 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[1]], collapse= "+")))
  }else{  form1 <-  as.formula(paste(name.response,' ~bols(Intercept,intercept=F)'))}
  
  if(any(RiskRedOver$parameter == parameter[2])){
    form2 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[2]], collapse= "+"))) 
  }else{  form2 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[3])){
    form3 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[3]], collapse= "+"))) 
  }else{  form3 <-  as.formula(paste(name.response, '~1'))}
  
  if(any(RiskRedOver$parameter == parameter[4])){
    form4 <- as.formula(paste(name.response, " ~ ", paste(RiskRedOver$VarName[RiskRedOver$parameter == parameter[4]], collapse= "+"))) 
  }else{  form4 <-  as.formula(paste(name.response, '~1,'))}
  
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
      model_after = gamboostLSS(formula, data = data, families = fam, method = 'noncyclic', weights = model.weights(object), control = boost_control(trace = T, mstop =  mstop(object), nu = list(as.numeric(object[[1]]$control$nu),as.numeric(object[[2]]$control$nu),as.numeric(object[[3]]$control$nu),as.numeric(object[[4]]$control$nu),as.numeric(object[[5]]$control$nu))))
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


data = attr(gauss_tanh_boost_fit_final,'data')
mod_desel <- DeselectBoostLSS_5(gauss_tanh_boost_fit_final, fam = Gauss_Cop(marg1 = "LOGLOG",
                                                                            marg2 = "LOGLOG",
                                                                            sigma1 = sd(data[model.weights(gauss_tanh_boost_fit_final)==1,]$GEBLAENGEK),
                                                                            sigma2 = sd(data[model.weights(gauss_tanh_boost_fit_final)==1,]$GEWICHT_sca),
                                                                            stabilization = "L2"),tau = 0.01)


Coefficients_desel <- mod_desel$coef


data = attr(gauss_tanh_boost_fit_final,'data')
test = data[model.weights(gauss_tanh_boost_fit_final)==0,]
test$Intercept <- 1

#Prediction error
pred_mu1_desel    <- predict(mod_desel$mu1, newdata = test, type = 'response')
pred_sigma1_desel <- predict(mod_desel$sigma1, newdata = test, type = 'response')
pred_mu2_desel    <- predict(mod_desel$mu2, newdata = test, type = 'response')
pred_sigma2_desel <- predict(mod_desel$sigma2, newdata = test, type = 'response')
pred_rho_desel    <- predict(mod_desel$rho, newdata = test, type = 'response')


LogLik_desel <- sum(loss(y = test[,c("GEBLAENGEK", "GEWICHT_sca")], mu1 = pred_mu1_desel, mu2 = pred_mu2_desel,
                         sigma1 = pred_sigma1_desel, sigma2 = pred_sigma2_desel, rho = pred_rho_desel), na.rm = T)





########################################
############ Stabsel  ##################
gauss_tanh_boost_fit <- gamboostLSS(eq_list,
                                    data = data_subset,
                                    weights = weights,
                                    families =Gauss_Cop(marg1 = "LOGLOG",
                                                        marg2 = "LOGLOG",
                                                        sigma1 = sd(data_subset[weights==1,]$GEBLAENGEK),
                                                        sigma2 = sd(data_subset[weights==1,]$GEWICHT_sca),
                                                        stabilization = "L2"),
                                    method = "noncyclic",
                                    control = boost_control(mstop = 8000,
                                                            nu = 0.01,
                                                            trace = T))



q = 20
pfer = 5

s1 <- stabsel(gauss_tanh_boost_fit, q = q, PFER = pfer, assumption = 'none')
sel = names(s1$selected)  


if(any(sel %in% (grep("(Intercept)",sel , value = T)))){
  sel <- sel[-which(sel %in% (grep("(Intercept)", sel, value = T)))]
}


data_subset$Intercept <- 1
if(length(name.response) > 1) name.response <- paste("cbind(",paste(name.response, collapse = ","),")")

name.response <- "cbind(GEBLAENGEK, GEWICHT_sca)"
if(any(grepl('.mu1', sel))){
  sel_mu1 <- gsub(".mu1","\\1",sel[grep('mu1',sel)])
  form_mu1 <- as.formula(paste(name.response,"~",paste(sel_mu1,collapse = '+'), sep= ''))
}else{
  sel_mu1 = 0
  form_mu1 <- as.formula(paste(name.response,"~", 'bols(Intercept, intercept=F)'))
}

if(any(grepl('mu2', sel))){
  sel_mu2 <- gsub(".mu2","\\1",sel[grep('mu2', sel)])
  form_mu2 <- as.formula(paste(name.response,"~",paste(sel_mu2,collapse = '+'),sep = ''))
}else{
  sel_mu2 = 0
  form_mu2 <- as.formula(paste(name.response,"~", 'bols(Intercept, intercept=F)'))
}

if(any(grepl('sigma1', sel))){
  sel_sigma1 <- gsub(".sigma1","\\1",sel[grep('sigma1', sel)])
  form_sigma1 <- as.formula(paste(name.response,"~",paste(sel_sigma1,collapse = '+'), sep = ''))
  
}else{
  sel_sigma1 = 0
  form_sigma1 <- as.formula(paste(name.response,"~", 'bols(Intercept, intercept=F)'))
}

if(any(grepl('sigma2', sel))){
  sel_sigma2 <- gsub(".sigma2","\\1",sel[grep('sigma2', sel)])
  form_sigma2 <- as.formula(paste(name.response,"~",paste(sel_sigma2,collapse = '+'),sep=''))
}else{
  sel_sigma2 = 0
  form_sigma2 <- as.formula(paste(name.response,"~", 'bols(Intercept, intercept=F)'))
}

if(any(grepl('rho', sel))){
  sel_rho <- gsub(".rho","\\1",sel[grep('rho', sel)])
  form_rho <- as.formula(paste(name.response,"~",paste(sel_rho,collapse = '+'),sep = ''))
}else{
  sel_rho = 0
  form_rho <- as.formula(paste(name.response,"~", 'bols(Intercept, intercept=F)'))
}

# abfrage, ob 
form <- list(mu1 = form_mu1,
             mu2 = form_mu2,
             sigma1 = form_sigma1,
             sigma2 = form_sigma2,
             rho = form_rho)



gauss_tanh_boost_fit <- gamboostLSS(form,
                                    data = data_subset,
                                    weights = weights,
                                    families =Gauss_Cop(marg1 = "LOGLOG",
                                                        marg2 = "LOGLOG",
                                                        sigma1 = sd(data_subset[weights==1,]$GEBLAENGEK),
                                                        sigma2 = sd(data_subset[weights==1,]$GEWICHT_sca),
                                                        stabilization = "L2"),
                                    method = "noncyclic",
                                    control = boost_control(mstop = 8000,
                                                            nu = 0.01,
                                                            trace = T))

set.seed(2)

folds <- cv(model.weights(gauss_tanh_boost_fit), type = "kfold", B = 10)
head(folds)
colSums(folds)

CV_mod <- cvrisk(gauss_tanh_boost_fit, folds = folds)

mstop(CV_mod)
plot(CV_mod)

# run final model on optimal mstop
gauss_tanh_boost_fit_final <- gamboostLSS(eq_list,
                                          data = data_subset,
                                          families = Gauss_Cop(marg1 = "LOGLOG",
                                                               marg2 = "LOGLOG", 
                                                               sigma1 = sd(data_subset_train$GEBLAENGEK),
                                                               sigma2 = sd(data_subset_train$GEWICHT_sca),
                                                               
                                                               stabilization = "L2"),
                                          weights = weights,
                                          method = "noncyclic",
                                          control = boost_control(mstop = mstop(CV_mod),
                                                                  nu = 0.01,
                                                                  trace = T))





test = data[model.weights(gauss_tanh_boost_fit_final)==0,]
test$Intercept <- 1

#Prediction error
pred_mu1_stabsel    <- predict(mod_desel$mu1, newdata = test, type = 'response')
pred_sigma1_stabsel <- predict(mod_desel$sigma1, newdata = test, type = 'response')
pred_mu2_stabsel    <- predict(mod_desel$mu2, newdata = test, type = 'response')
pred_sigma2_stabsel <- predict(mod_desel$sigma2, newdata = test, type = 'response')
pred_rho_stabsel    <- predict(mod_desel$rho, newdata = test, type = 'response')


LogLik_desel <- sum(loss(y = test[,c("GEBLAENGEK", "GEWICHT_sca")], mu1 = pred_mu1_stabsel, mu2 = pred_mu2_stabsel,
                         sigma1 = pred_sigma1_stabsel, sigma2 = pred_sigma2_stabsel, rho = pred_rho_stabsel), na.rm = T)



########################################
############ Probing  ##################

set.seed(1)

weights <- sample(c(rep(1, 4103), rep(0, 2000)))
data_subset_train <- data_subset[weights == 1, ]

data_subset %>% select(!c(GEBLAENGEK, GEWICHT_sca))
probes <- as.data.frame(sapply(data_subset %>% select(!c(GEBLAENGEK, GEWICHT_sca)), sample)) # shuffling all variables (leaving out response which is at position 1 and 2)
names(probes) <- paste(names(data_subset %>% select(!c(GEBLAENGEK, GEWICHT_sca))), "probe", sep = "_")

probes.train <- cbind(data_subset,probes)

### formula
names(data_subset)

mu1_eq <- cbind(GEBLAENGEK, GEWICHT_sca) ~ bols(GDiabetis) + 
  bols(GESCHLECHT) + 
  bbs(ABDOMEN.SAGITTALDURCHMESSER) + 
  bbs(ABDOMEN.TRANSVERSALDURCHMESSE) + 
  bbs(ABDOMENUMFANG) + 
  bbs(BIPARIETALER.DURCHMESSER) + 
  bbs(FRONTOOKZIPITALER.DURCHMESSER) + 
  bbs(KOPFUMFANG) + 
  bbs(FEMURLAENGE) +
  bbs(ABSAG_ABTRAN) +
  bbs(ABSAG_ABUMF) +
  bbs(ABSAG_BIPADU) +
  bbs(ABSAG_FRODU) +
  bbs(ABSAG_KOPUM) +
  bbs(ABSAG_FEMLA) +
  bbs(ABTRAN_ABDUM) +
  bbs(ABTRAN_BIPADU) + 
  bbs(ABTRAN_FRONDU) +
  bbs(ABTRAN_KOPUM) +
  bbs(ABTRAN_FEMLA) + 
  bbs(ABDUM_BIPADU) +
  bbs(ABDUM_FRONDU) +
  bbs(ABD_KOPUM) +
  bbs(ABDUM_FEMLA) +
  bbs(BIPADU_FRODU) +
  bbs(BIPADU_KOPUM) +
  bbs(BIPADU_FEMLA) +
  bbs(FRONDU_KOPUM) +
  bbs(FRONDU_FEMLA) +
  bbs(KOPUM_FEMLA) +
  bbs(ENDGEWICHT) + 
  bbs(MUTTERGROESSE) + 
  bbs(bmi_mutter) + 
  bbs(GRAVIDA) +
  bbs(PARA) +
  bbs(GestAge_Days_US) + ###
  bols(GDiabetis_probe) + 
  bols(GESCHLECHT_probe) + 
  bbs(ABDOMEN.SAGITTALDURCHMESSER_probe) + 
  bbs(ABDOMEN.TRANSVERSALDURCHMESSE_probe) + 
  bbs(ABDOMENUMFANG_probe) + 
  bbs(BIPARIETALER.DURCHMESSER_probe) + 
  bbs(FRONTOOKZIPITALER.DURCHMESSER_probe) + 
  bbs(KOPFUMFANG_probe) + 
  bbs(FEMURLAENGE_probe) +
  bbs(ABSAG_ABTRAN_probe) +
  bbs(ABSAG_ABUMF_probe) +
  bbs(ABSAG_BIPADU_probe) +
  bbs(ABSAG_FRODU_probe) +
  bbs(ABSAG_KOPUM_probe) +
  bbs(ABSAG_FEMLA_probe) +
  bbs(ABTRAN_ABDUM_probe) +
  bbs(ABTRAN_BIPADU_probe) + 
  bbs(ABTRAN_FRONDU_probe) +
  bbs(ABTRAN_KOPUM_probe) +
  bbs(ABTRAN_FEMLA_probe) + 
  bbs(ABDUM_BIPADU_probe) +
  bbs(ABDUM_FRONDU_probe) +
  bbs(ABD_KOPUM_probe) +
  bbs(ABDUM_FEMLA_probe) +
  bbs(BIPADU_FRODU_probe) +
  bbs(BIPADU_KOPUM_probe) +
  bbs(BIPADU_FEMLA_probe) +
  bbs(FRONDU_KOPUM_probe) +
  bbs(FRONDU_FEMLA_probe) +
  bbs(KOPUM_FEMLA_probe) +
  bbs(ENDGEWICHT_probe) + 
  bbs(MUTTERGROESSE_probe) + 
  bbs(bmi_mutter_probe) + 
  bbs(GRAVIDA_probe) +
  bbs(PARA_probe) +
  bbs(GestAge_Days_US_probe) 

sigma1_eq <- mu2_eq <- sigma2_eq <- rho_eq <- mu1_eq

eq_list <- list(mu1 = mu1_eq,
                sigma1 = sigma1_eq,
                mu2 = mu2_eq,
                sigma2 = sigma2_eq,
                rho = rho_eq)



####################### Gauss: tanh link function ##############################

source("Copulas/Copula_Gaussian.R")

gauss_tanh_boost_fit <- gamboostLSS(eq_list,
                                    data = probes.train,
                                    weights = weights,
                                    families =  Gauss_Cop(marg1 = "LOGLOG",
                                                          marg2 = "LOGLOG",
                                                          sigma1 = sd(probes.train[weights==1,]$GEBLAENGEK),
                                                          sigma2 = sd(probes.train[weights==1,]$GEWICHT_sca),
                                                          stabilization = "L2"),
                                    method = "noncyclic",
                                    control = boost_control(mstop = 5000,
                                                            nu = 0.01,
                                                            trace = T))
p = 36
mstop_probes <- min(which(selected(gauss_tanh_boost_fit,merge = T) > p   )) - 1 

coef(gauss_tanh_boost_fit)
names(coef(gauss_tanh_boost_fit)$mu1)
names(coef(gauss_tanh_boost_fit)$sigma1)
names(coef(gauss_tanh_boost_fit)$mu2)
names(coef(gauss_tanh_boost_fit)$sigma2)
names(coef(gauss_tanh_boost_fit)$rho)

gauss_tanh_boost_fit[mstop_probes]
gauss_tanh_boost_fit_final <- gauss_tanh_boost_fit


#Prediction error
pred_mu1_prob    <- predict(gauss_tanh_boost_fit$mu1, newdata = test, type = 'response')
pred_sigma1_prob <- predict(gauss_tanh_boost_fit$sigma1, newdata = test, type = 'response')
pred_mu2_prob    <- predict(gauss_tanh_boost_fit$mu2, newdata = test, type = 'response')
pred_sigma2_prob <- predict(gauss_tanh_boost_fit$sigma2, newdata = test, type = 'response')
pred_rho_prob    <- predict(gauss_tanh_boost_fit$rho, newdata = test, type = 'response')

LogLik_prob <- sum(loss(y = test[,c("GEBLAENGEK", "GEWICHT_sca")], mu1 = pred_mu1_prob, mu2 = pred_mu2_prob,
                         sigma1 = pred_sigma1_prob, sigma2 = pred_sigma2_prob, rho = pred_rho_prob), na.rm = T)






