#--------------------------------------------------------------------------------
## Title:                       Main Script for a Binary Simulation Study
## Author:                      Ermioni Athanasiadi
## Last modification:           21.08.2025   
## History of modifications:
##                              10.07.2025
## a simulation study with 
## binary outcome and simple one-factor experimental design
## manipulated conditions:
##      -- sample size n
##      -- effect size (coefficients: or, interc)
#--------------------------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
source("functions.R")

### REMOVE THIS AND ADD THE FOLLOWING TO functions.R in each function: 
# add ./data in file path of where to save or retrieve files from
setwd("./data")

n_vec <- c(20)
or_vec <- c(1.1)
interc <- 1.1
miss_vec <- c(0.5)

impmethod_vec <- c("none", "logreg")
impmethod <- "none"

nreps <- 10
seed = 2025

generate_replications(n_vec, or_vec, interc, nreps, seed)

analyse_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete = TRUE)
setwd("./../results")
summarize_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed, is.complete = TRUE)

setwd("./../data")
generate_na_data(n_vec, or_vec, interc, nreps, miss_vec, seed)

impute_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed)

analyse_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete = FALSE)

setwd("./../results")
summarize_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed, is.complete = FALSE)


# setwd("./../results")