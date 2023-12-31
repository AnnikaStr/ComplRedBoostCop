### defs for sims ###

library(batchtools)
library(gamboostLSS)
library(stabs)

OVERWRITE = FALSE

#settings data generation
p = c(200)
#p = c(1000)
n = c(1000)
n_test = 1000
n_eval = 1500
n_train = 1000

# deselection
TAU_1 = c(0.001,0.01)

# setting stabsel
Q = c(10, 15, 20)

# settings cluster
MEMORY =  1024L*20L
NTASKS = 10
NCPUS = 10
WALLTIME = 600L*600L

REPLICATIONS = 100
