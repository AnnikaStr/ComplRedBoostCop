# Deselection for boosting distributional copula regression
Function for complexity reduction of boosting distributional copula regression based on the deselection of base-learner.

# Example
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



