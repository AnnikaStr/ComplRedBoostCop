# Complexity reduction and enhanced variable selection for boosting distributional copula regression

## Example
```r
library(copula)
library(gamboostLSS)
source("Copulas/Copula_Gaussian.R") # source Gaussian copula function

n_train = 1000
n_eval = 1500
p = 20

set.seed(1)

# generating data
data.gen <- function(FAM, Z, i){
  
  
  z1 <- Z[i, 1]
  z2 <- Z[i, 2]
  z3 <- Z[i, 3]
  z4 <- Z[i, 4]
  
  
  eta_mu1    <- -1 * z1 + 0.5 * z3
  eta_sigma1 <- - 0.7  - 0.7 * z3
  eta_mu2    <- - 0.5 - 0.7 * z1 + 0.3 * z2
  eta_sigma2 <- 2 + 0.5 * z2
  eta_theta  <- 1 + 1 * z4
  
  if( FAM == "clayton") theta.para <- exp(eta_theta) + 1e-07
  if( FAM == "gumbel") theta.para <- exp(eta_theta) + 1
  if( FAM %in% c ("normal" , "t")) theta.para <- tanh(eta_theta)         
  if( FAM %in% c ("clayton" , "gumbel") ) Cop <- archmCopula(family = FAM,
                                                             dim = 2, param = theta.para)
  else Cop <- ellipCopula(family = FAM, dim = 2, param = theta.para , df = 4)
  speclist1 <- list( meanlog = eta_mu1, sdlog = exp(eta_sigma1) )
  speclist2 <- list( mu = exp(eta_mu2), sigma = exp(eta_sigma2), nu = 1, tau = 1)
  spec <- mvdc(copula = Cop, c("lnorm" , "GB2") , list(speclist1 , speclist2) )
  c(rMvdc(1, spec), Z[i, ])
}


data_train = as.data.frame(matrix(0, ncol = p+2+1, nrow=(n_train+n_eval)))
names(data_train) = c("y1", "y2", paste0("x", 1:p),'Intercept')
X_train = matrix(NA, ncol = p, nrow = (n_train+n_eval))

for(i in 1:p){
  X_train[,i] = runif((n_train+n_eval), -1, 1)
}

X_train = cbind(X_train, 1)

weights_train = c(rep(1, times = n_train), rep(0, times = n_eval))

for(i in 1:(n_train + n_eval)){
  data_run = data.gen(FAM = "normal", Z = X_train, i = i)
  data_train[i, ] = data_run
}


# fit a boosting model
mod = glmboostLSS(cbind(y1,y2) ~.,
                  data = data_train,
                  control = boost_control(mstop = 5000,
                                          nu = 0.01,
                                          risk = "oobag",
                                          trace=T),
                  weights = weights_train,
                  method = 'noncyclic',
                  families = Gauss_Cop(marg1 = "LOGNO", marg2 = "LOGLOG"))

MSTOP = which.min(risk(mod,merge = T))

oobag.risk = risk(mod,merge = T)

rm(mod) # removed the fist fitted model
data_train = data_train[weights_train == 1, ]

mod = glmboostLSS(cbind(y1,y2) ~.,
                  data = data_train,
                  control = boost_control(mstop = MSTOP,
                                          nu = 0.01),
                  method = 'noncyclic',
                  families = Gauss_Cop(marg1 = "LOGNO",marg2 = "LOGLOG"))



## Deselection with a threshold value of 1%
tau = 0.01
mod_desel =  DeselectBoost(mod, tau =.01, fam = Gauss_Cop(marg1 = "LOGNO",
                                                          marg2 = "LOGLOG"))


## Stability selection
set.seed(seed)
q = 20
pfer = 5

train1 <- train[weights_train == 1,]
  
mod = gamboostLSS(cbind(y1,y2) ~., data = train1,
                    control = boost_control(mstop = 5000, nu = 0.01), method = 'noncyclic', families = Gauss_Cop(marg1 = "LOGNO",
                                                                                                                marg2 = "LOGLOG"))


s1 <- stabsel(mod, q = q, PFER = pfer, assumption = 'none')
sel = names(s1$selected)  

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
 
    form <- list(mu1 = form_mu1,
                 mu2 = form_mu2,
                 sigma1 = form_sigma1,
                 sigma2 = form_sigma2,
                 rho = form_rho)
    
mod = gamboostLSS(form, data = train,
                      control = boost_control(mstop = 5000, nu = 0.01,risk = 'oobag'), weights = weights_train, method = 'noncyclic', families = Gauss_Cop(marg1 = "LOGNO",
                                                                                                                                                           marg2 = "LOGLOG")
MSTOP <- which.min(risk(mod,merge = T))
oobag.risk <- risk(mod,merge = T)

mod = gamboostLSS(form, data = train1, control = boost_control(mstop = MSTOP, nu = 0.01), method = 'noncyclic', families = Gauss_Cop(marg1 = "LOGNO",marg2 = "LOGLOG"))
    
## Probing
train <- train[weights_train == 1,]

set.seed(seed)

# shuffling all variables (leaving out response which is at position 1 and 2)
probes.train <- as.data.frame(sapply(train[,-c(1,2)], sample)) 
names(probes.train) <- paste(names(train)[-c(1,2)], "probe", sep = "_")
probes.train <- cbind(train,probes.train)

mod = gamboostLSS(cbind(y1,y2) ~., data = probes.train, control = boost_control(mstop = 2000, nu=0.01),
                    method = 'noncyclic', families = Gauss_Cop(marg1 = "LOGNO",
                                                               marg2 = "LOGLOG"))

if(!any(selected(mod,merge = T) > p)) mod[5000]
if(!any(selected(mod,merge = T) > p)) mod[10000]
mstop_probes <- min(which(selected(mod,merge = T) > p   )) - 1 
mod[mstop_probes]

