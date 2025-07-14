#--------------------------------------------------------------------------------
## Title:                       Main Script for a Binary Simulation Study
## Author:                      Ermioni Athanasiadi
## Last modification:           10.07.2025
## a simulation study with 
## binary outcome and simple one-factor experimental design
## manipulated conditions:
##      -- sample size n
##      -- effect size (coefficients: or, interc)
#--------------------------------------------------------------------------------

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
source("functions.R")

setwd("./data")

n_vec <- c(50, 100)
or_vec <- c(1.1, 2)
interc <- 1.1
miss_vec <- c(0.5, 0.25)
impmethod <- "pmm"

nreps <- 50
#single_sim <- generate_df(50, or, interc)
seed = 2025

generate_replications(n_vec, or_vec, interc, nreps, seed)

generate_na_data(n_vec, or_vec, interc, nreps, miss_vec, seed)

impute_data(n_vec, or_vec, interc, nreps, missr, impmethod, seed)

analyse_data(n_vec, or_vec, interc, nreps, missr, impmethod, seed)
