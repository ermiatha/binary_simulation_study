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

n_vec <- c(20, 30, 50, 80, 100, 130, 150)
#n_vec <- c(20)
or_vec <- c(1.1, 1.6, 2)
interc <- 1.1
miss_vec <- c(0.1, 0.25, 0.5)

impmethod_vec <- c("pmm", "logreg")
nreps <- 500
seed = 2025


### PMM  ###############
impmethod <- "pmm"
impmethod_vec <- c("pmm")

generate_replications(n_vec, or_vec, interc, nreps, seed)

# analyse_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete = TRUE)
# -> Work in progress: rewrite as function analyze_complete_data()
# setwd("./../results")
# summarize_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed, is.complete = TRUE)

setwd("./../data")
generate_na_data(n_vec, or_vec, interc, nreps, miss_vec, seed)

impute_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed)

analyse_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete = FALSE)

setwd("./../results")
summarize_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed, is.complete = FALSE)


# setwd("./../results")

setwd("./../data")
### LOGREG  ###############
impmethod <- "logreg"
impmethod_vec <- c("logreg")

generate_replications(n_vec, or_vec, interc, nreps, seed)

# analyse_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete = TRUE)
# -> Work in progress: rewrite as function analyze_complete_data()
# setwd("./../results")
# summarize_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed, is.complete = TRUE)

setwd("./../data")
generate_na_data(n_vec, or_vec, interc, nreps, miss_vec, seed)

impute_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed)

analyse_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete = FALSE)

setwd("./../results")
summarize_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed, is.complete = FALSE)


##

setwd("./../data")
### BLIMP  ###############
impmethod <- "blimp"
impmethod_vec <- c("blimp")

generate_replications(n_vec, or_vec, interc, nreps, seed)

# analyse_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete = TRUE)
# -> Work in progress: rewrite as function analyze_complete_data()
# setwd("./../results")
# summarize_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed, is.complete = TRUE)

setwd("./../data")
generate_na_data(n_vec, or_vec, interc, nreps, miss_vec, seed)

impute_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed)

analyse_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete = FALSE)

setwd("./../results")
summarize_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed, is.complete = FALSE)


setwd("./../data")
### Listwise Deletion  ###############
impmethod <- "ld"
impmethod_vec <- c("ld")

analyse_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete = FALSE)

setwd("./../results")
summarize_data(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed, is.complete = FALSE)



##############################################
# visualize pmm results

impmethod <- "pmm"
visualize_bias(impmethod, show_title = FALSE)
visualize_power(impmethod, show_title = FALSE)
visualize_coverage(impmethod, show_title = FALSE)
visualize_convergence(impmethod, show_title = FALSE)
visualize_runtime(impmethod, show_title = FALSE)
visualize_runtime_by_impmethod(show_title = FALSE)


# visualize logreg results
impmethod <- "logreg"
visualize_bias(impmethod, show_title = FALSE)
visualize_power(impmethod, show_title = FALSE)
visualize_coverage(impmethod, show_title = FALSE)
visualize_convergence(impmethod, show_title = FALSE)
visualize_runtime(impmethod, show_title = FALSE)
visualize_runtime_by_impmethod(show_title = FALSE)


# visualize LD results
impmethod <- "ld"
visualize_bias(impmethod, show_title = FALSE)
visualize_power(impmethod, show_title = FALSE)
visualize_coverage(impmethod, show_title = FALSE)
visualize_convergence(impmethod, show_title = FALSE)
visualize_runtime(impmethod, show_title = FALSE)
visualize_runtime_by_impmethod(impmethods = c("pmm", "logreg"), show_title = FALSE)






