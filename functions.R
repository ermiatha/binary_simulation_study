
#--Binary simulation study-------------------------------------------------------
## Title:                       Functions for a Binary Simulation Study
## Author:                      Ermioni Athanasiadi
## Last modification:           10.07.2025
## binary outcome and simple one-factor experimental design
## manipulated conditions:
## manipulated conditions:
##      -- sample size n
##      -- effect size (coefficients: or, interc)


#--Setup ------------------------------------------------------------------------
# this does not run from the console  # set working directory to where file is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(dplyr)
library(ggplot2)
library(mice)
library(future.apply)
library(rblimp)

#________________________________________________________________________________
# function to check valid user input for sample size
check_even <- function(n) {
    if (any(n %% 2 != 0)) {
        stop("Please select even number as sample size")
    } 
}


#-Part I: Data generation-------------------------------------------------------
# Function to generate Data for single dataset
generate_df <- function(n, or, interc) {
    
    # check for valid input of n
    check_even(n)
    
    # define variables
    group <- factor(rep(0:1, each = n/2)) # binary treatment
    id <- seq(1:n) # id
    # design matrix
    X <- model.matrix(~ group)  # contains only the fixed effect
    # define coefficients
    beta <- numeric(ncol(X))
    beta[1] <- log(interc) # intercept: logg odds for add in control group # gives prob = 0.5 
    beta[2] <- log(or) # diff of treatment group to control group
    # convert to probability
    true_prob <- exp(beta[2]) / (1 + exp(beta[2]))
    # linear outcome then binary
    eta_i <- (X %*% beta)
    y_i <- rbinom(nrow(X), size = 1, prob = binomial()$linkinv(eta_i))
    list(id = id, group = group, y_i = y_i, true_prob = true_prob, X = X, beta = beta)
    
    # Generate data and create dataframe
    dat <- data.frame(id = id, group = group, y_i = as.factor(y_i))
    
    # Track amount of zeros in outcome
    perc_zeros <- length(which(dat$y_i == 0)) / length(dat$y_i)
    
    return(dat)
}

#________________________________________________________________________________
# Function to create replications of generated datasets
generate_replications <- function(n_vec, or_vec, interc, nreps, seed) {

    
    # set seed
    set.seed(seed)
    
    # Double check: valid number of n?
    # sapply(n_vec, check_even)
    
    
    for (j in 1:length(n_vec)) {
        n <- n_vec[j]
        
        for (i in 1:length(or_vec)) {
            or <- or_vec[i]
            # set variable name for use in global environment
            vname <- paste('list_data',n,'_',or, sep='')
            print(paste("current list:", vname))
            time.parallel <- system.time(
                list_sims <- future_lapply(
                    1:nreps, function(i) generate_df(
                        n = n, 
                        or = or, 
                        interc = interc
                    ), 
                    future.seed = TRUE))
            # assign each list to a different name
            assign(vname, list_sims)
            fname <- paste('list_data',n,'_',or,'.RData', sep='')
            save(list_sims, file = fname, compress = F)
            
        }
    }
}

# export execution time ###
# time.parallel

#--Part II: NA Generation--------------------------------------------------------
# Function to create missing values given a dataframe and a missingness rate
generate_na <- function(dat, missr) {
    
    # Track amount of zeros in outcome
    perc_zeros <- length(which(dat$y_i == 0)) / length(dat$y_i)
    
    # Create MAR missingness
    missr <- missr * 2  # e.g. 40% of N/2 sample --> 20% of complete sample
    group0_indices <- which(dat$group == 0)
    group1_indices <- which(dat$group == 1)

    
    # Create MAR missingness: allows a total missr of up to 60%
    miss_indeces0 <- sample(group0_indices, size = floor(length(group0_indices) * missr*0.8))
    miss_indeces1 <- sample(group1_indices, size = floor(length(group0_indices) * missr*0.2))
    dat$y_i[miss_indeces0] <- NA
    dat$y_i[miss_indeces1] <- NA
    
    return(dat)
}

#________________________________________________________________________________
# Function to apply na generating function to all replications and conditions

generate_na_data <- function(n_vec, or_vec, interc, nreps, miss_vec, seed) {
    for (j in 1:length(n_vec)) {
        n <- n_vec[j]
        
        for (i in 1:length(or_vec)) {
            or <- or_vec[i]
            fname <- paste('list_data',n,'_',or,'.RData', sep='')
            
            # get object of loaded file
            list_sims <- get(load(fname) )  
            
            # create different missingness levels 
            for (m in 1:length(miss_vec)) {
                missr <- miss_vec[m]
                
                ## introduce missings in each replication
                list_sims <- future_lapply(list_sims, function(i) generate_na(i, missr = missr), future.seed = TRUE)
                
                fname <- paste('list_data_na',n,'_',or,'_',missr,'.RData', sep='')
                vname <- paste('list_data_na',n,'_',or,'_',missr, sep='')
                
                assign(vname, list_sims)
                
                # save incomplete datasets as files
                print(paste("Saving file:", vname))
                save(list_sims, file = fname, compress = F)
            }
        }
    } 
}

#--Part III: Data Imputation / Handling ----------------------------------------
# Function to impute missing values with blimp given a dataframe
impute_blimp <- function(dat) {
    
    dat$y_i <- as.integer(dat$y_i)
    dat$y_i <- ifelse(dat$y_i == 1, 0, 1)
    
    # truemodel <- glm(y_i ~ group, data = dat, family = binomial())
    # summary(truemodel)
    
    imp_blimp <- rblimp(
        data = dat,
        ordinal = 'y_i group',
        fixed = 'group',
        model = 'logit(y_i) ~ group',
        seed = 2025,
        burn = 1000,
        iter = 10000,
        nimps = 5
    )
    
    # output(imp_blimp)
    imp <- as.mitml(imp_blimp)
    

    return(imp)
}

#________________________________________________________________________________
# Function to impute missing values with mice given a dataframe
impute_mice <- function(dat, method) {
    dat <- dat %>% select(group, y_i)
    meth <- make.method(dat)
    meth["y_i"] <- method
    pred <- make.predictorMatrix(dat)
    
    # tryCatch to handle any errors and return NULL if imputation fails
    imp <- tryCatch({
        suppressWarnings({
            mice(dat, method = meth, predictorMatrix = pred, m = 5, print = F)
        })
    }, error = function(e) {
        message("Error in imputation: ", e$message)  # Print the error message
        return(NULL)  # Return NULL to indicate failure
    }) 
    
    return(imp)
}

#________________________________________________________________________________
# Function to impute missing values given a dataframe and a method
impute_method <- function(dat, method) {
    
    if (method == "blimp") {
        imp <- impute_blimp(dat)
    }
    
    else if (method == "pmm" | method == "logreg") {
        imp <- impute_mice(dat, method)
    }
    
    return(imp)
}

#________________________________________________________________________________
# Function to perform listwise deletion

ld_method <- function(dat) {
    dat <- na.omit(dat)
    return(dat)
}


#________________________________________________________________________________
# Function to generate imputed data for all replications and conditions
impute_data <- function(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed) {
    for (j in 1:length(n_vec)) {
        n <- n_vec[j]
        
        for (i in 1:length(or_vec)) {
            or <- or_vec[i]
            
            
            # create different missingness levels 
            for (m in 1:length(miss_vec)) {
                
                # load dataframes
                missr <- miss_vec[m]
                fname <- paste('list_data_na',n,'_',or,'_',missr,'.RData', sep='')
                list_sims <- get(load(fname) )  # get object of loaded file
                # print(fname)

                # if valid impmethod, generate imputated datasets for all replications
                if (impmethod != "none") {
                    # generate imputated datasets for all replications

                    time.parallel <- system.time(
                        list_sims <- future_lapply(
                            list_sims, function(i) impute_method(
                                i, 
                                method = impmethod),
                            future.seed = TRUE))
                    
                    # print execution time for imputation
                    print(time.parallel)    
                }

                
                ## set filepath for result folder
                filepath = "./../../binary_simulation_study/results/"
                
                # save convergence information as df to file #
                ##
                ## 
                ##
                ##
                # future idea:
                # only variable to be retrieved in this step, all others later
                # --> rename for clarity ("convergence"), save as vectors?
                ##
                fname <- paste('df_perf_results',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                df_perf_results <- performance_results(list_sims, nreps)
                print(paste("Saving file:", fname))
                save(df_perf_results, file = paste0(filepath, fname), compress = F)
                
                ## set file names and save raw data to files
                fname <- paste('list_data_imp',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                vname <- paste('list_data_imp',n,'_',or,'_',missr,impmethod, sep='')
                assign(vname, list_sims)
                print(paste("Saving file:", vname))
                save(list_sims, file = fname, compress = F)
                
                
                # save time information as vector to results / data folder HERE
                ##
                ##
                ###
            }
        }
    }
}

#--Part IV: Data Analysis--------------------------------------------------------
# Function to get performance_results:
performance_results <- function(list_sims, nreps) {
    # Function to save results about performance:
    # ----- convergence
    # ----- runtime
    
    # set up dataframe for all nreps replications
    df <- data.frame(matrix(nrow = nreps, ncol = 1, dimnames = list(c(), c("convergence"))))  #  
    # "runtime"  not needed here, only global runtime over all nreps available due to parallelization
    
    # create dummies for convergence performance
    list_convergence <- cbind(lapply(list_sims, typeof))
    ## save convergence information to list: 1 indicates success, 0 (NULL) failure
    df$convergence <- ifelse(unlist(list_convergence) == "list", 1, 0)
    
    return(df)
}

#________________________________________________________________________________
# function to estimate effects for single dataframe given a list of imputed datasets
estimate_effects_imp <- function(imp) {
    
    
    # initiate empty results df
    model_coefs <- data.frame("estimate" = NA, "std.error" = NA, "statistic" = NA, 
                              "df" = NA, "p.value" = NA, "conf.low" = NA, "conf.high" = NA)
    
        # only estimate pooled results if replication is successful
        if (!is.null(imp)) {
            imp_models <- with(data = imp, glm(y_i ~ group, family = "binomial"))
            # save only result for coefficient of group effect, ignore intercept estimate
            model_coefs <- summary(mice::pool(imp_models), conf.int = TRUE)[-1, -1]
            # remove double columns for lower and upper CI
            model_coefs <- model_coefs[-c(6,7)]

        }

    return(model_coefs)
}

#________________________________________________________________________________
# function to estimate effects for single dataframe given a complete or incomplete dataset
estimate_effects <- function(df) {
    # initiate empty results df
    model_coefs <- data.frame("estimate" = NA, "std.error" = NA, "statistic" = NA, 
                              "df" = NA, "p.value" = NA, "conf.low" = NA, "conf.high" = NA)
    
    # estimate effects with incomplete dataset -> handles missings with LD
    model_summ <- summary(glm(y_i ~ group, family = "binomial", data = df))$coefficients[-1, ]
    
    # extract relevant parameters for output
    # MISSING INFORMATION FOR LOWER AND UPPER CI
    model_coefs$estimate <- model_summ[1]
    model_coefs$std.error <- NULL
    model_coefs[, c(1, 2, 4)] <- c(model_summ[1], model_summ[2], model_summ[4])
    
    
    return(model_coefs)
}

#________________________________________________________________________________
# Function to calculate performance-related outcomes, given a dataframe
add_perf_results <- function(df, or)  { #what other parameters? 
    # given: a dataframe containing all rows with successful replications
    
    true_or <- or
    
    # boolean: significant p-value y=1 or n=0
    df$signif <- ifelse(df$p.value < .05, 1, 0)
    
    df$coverage_dummy <- ifelse((true_or >= df$conf.low & true_or <= df$conf.high), 1, 0)
    
    # Bias
    # abs_bias for or
    # rel_bias for or
    
    # abs_bias for se
    # rel_bias for se
    
    ## percent zeros
    
    ## runtimes ##
    
    return(df)
}

#________________________________________________________________________________
# for future restructuring: migrate analysis for complete data into extra function!
#
# Function to analyze complete data as benchmark for performance of imputation methods
# analyse_complete_data <- function(n_vec, or_vec, interc, nreps, miss_vec, seed) {
#     
# }
# 


#________________________________________________________________________________
# Function to analyse data and export dfs with results for all replications and conditions 
analyse_data <- function(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete) {
    
    for (j in 1:length(n_vec)) {
        n <- n_vec[j]
        
        for (i in 1:length(or_vec)) {
            or <- or_vec[i]
            
            # different missingness levels 
            for (m in 1:length(miss_vec)) {
                print(paste("Currently at n, or, missingness:", n_vec[j], or_vec[i], miss_vec[m]))
                
                # load imputed dataframes
                missr <- miss_vec[m]
                fname <- paste('list_data_imp',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                list_sims <- get(load(fname) )  # get object of loaded file
                
                # specify filepath for retrieving and saving results
                filepath_results = "./../../binary_simulation_study/results/"
                # retrieve df with performance results
                fname_perf <- paste('df_perf_results',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                
                df_perf_results <- get(load(paste0(filepath_results, fname_perf)))
                
                # if no impmethod specified or data is complete: simple effect estimation
                # Note: this performs listwise deletion if data is incomplete
                if (impmethod == "none" | is.complete == TRUE) {
                    time.parallel <- system.time(
                    estimated_effects <- future_lapply(list_sims, function(i) estimate_effects(i), future.seed = TRUE)
                )
                }
                
                # if impmethod specified: estimate effects with imp object
                else {
                    # estimate effects for all replications (NA rows for non-converged replications)
                    time.parallel <- system.time(
                        estimated_effects <- future_lapply(list_sims, function(i) estimate_effects_imp(i), future.seed = TRUE)
                    )
                }
                
                ## user info:
                num_successful_sims = length(which(list_sims != "NULL"))
                perc_successful_sims = (num_successful_sims / length(list_sims) ) * 100
                print(paste("Effects were estimated based on successful replications: 
                            Number of successful replications:", 
                            num_successful_sims, 
                            "Percentage % ", 
                            perc_successful_sims))

                
                # print time for extracting results
                print(time.parallel)

                # save results of all replications as dataframes
                df_results_all <- do.call(rbind, lapply(estimated_effects, function(x) x))
                # add performance results to df
                df_results_all <- cbind(df_results_all, df_perf_results)
                
                # give different filename to complete dataframe
                ##
                # future idea: migrate data analysis for complete data into extra 
                # function for clarity and structure
                ##
                ##
                if (is.complete == TRUE) {
                    print("correctly checking that we have a complete dataset here")
                    complete <- "complete"
                }
                else {
                    complete <- ""
                }
                
                # set file names
                fname <- paste(complete,'df_results_all',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                vname <- paste(complete,'df_results_all',n,'_',or,'_',missr,impmethod, sep='')
                assign(vname, df_results_all)
                # save complete result dfs as files
                print(paste("Saving file:", vname))
                save(df_results_all, file = paste0(filepath_results, fname), compress = F)
                
                # now look at subset of successful replications (convergence)
                # extract relevant rows 
                df_results_succ <- df_results_all %>%
                    filter(convergence == 1)
                
                # add further variables for performance-related outcomes #
                df_results_succ <- add_perf_results(df_results_succ, or)
                
                # set file names
                fname <- paste(complete,'df_results_succ',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                vname <- paste(complete,'df_results_succ',n,'_',or,'_',missr,impmethod, sep='')
                assign(vname, df_results_succ)
                # save successful result dfs as files
                print(paste("Saving file:", vname))
                save(df_results_succ, file = paste0(filepath_results, fname), compress = F)

                
            }
        }
    }
}

#--Part VI: Summary of Results -------------------------------------------------
# Function to summarize the results of successful replications of a simulation study, given a df
get_summary_stats <- function(df, n, or, missr, nreps) {
    # input: a dataframe with results for all successful replications, other pars
    
    mean_or_hat   <- mean(exp(df[, 1]))
    or_hat        <- exp(df[, 1])
    mean_prob_hat <- mean(exp(df[, 1])) / (1 + mean(exp(df[, 1])))
    power         <- sum(df$signif) / nrow(df)
    true_or       <- or
    true_prob     <- or / (1 + or)   # check this!!
    coverage      <- sum(df$coverage_dummy) / nrow(df)
    perc_convergence <- sum(df$convergence) / nreps * 100
    
    # absolute and relative bias for or
    abs_bias_or   <- mean_or_hat - true_or
    rel_bias_or   <- mean(abs(true_or - or_hat) / true_or * 100)
    
    # absolute and relative bias for standard error
    mean_se_hat   <- mean(df[, 2])
    se_hat        <- df[, 2]
    # true_se has to be calculated earlier based on complete sample and then passed
    # as argument to this function
    true_se       <- 1.1111111111111111111111111
    abs_bias_se      <- mean_se_hat - true_se
    rel_bias_se      <- mean(abs(true_se - se_hat) / true_se * 100)
    
    # perc_zeros
    
    df_summ <- data.frame(
                            mean_or_hat, mean_prob_hat, 
                            true_or, true_prob,
                            mean_se_hat,
                            power, coverage,
                            perc_convergence
                          )

    return(df_summ)
    
}



#________________________________________________________________________________
# Function to summarize and export results as dfs for all simulation conditions
summarize_data <- function(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed, is.complete) {

    # extend this later to run for all impmethods ? ###
    impmethod <- impmethod_vec[1]
    ##
    ##
    ##
    
    if (is.complete == TRUE) {
        complete <- "complete"
    }
    else {
        complete <- ""
    }
    
    # specify filepath for retrieving and saving results
    filepath_results = "./../results/"
    fname <- paste(complete,'df_results_succ',n_vec[1],'_',or_vec[1],'_',miss_vec[1],impmethod,'.RData', sep='')
    df_results_succ <- get(load(fname) )  # get object of loaded file
    
    # initialize empty result dataframe for all conditions (except impmethod, see above)
    df_summ_all_cond <- expand.grid(impmethod = impmethod, 
                                    miss = miss_vec, 
                                    or = or_vec, 
                                    n = n_vec)[, c("n", "or", "miss", "impmethod")]
    # grab all result variables 
    result_vars <- names(get_summary_stats(df = df_results_succ, n = n_vec[1], or_vec[1], miss_vec[1], nreps))
    df_summ_all_cond[, result_vars] <- NA
    
    # get results for each condition
    for (j in 1:length(n_vec)) {
        n <- n_vec[j]
        
        for (i in 1:length(or_vec)) {
            or <- or_vec[i]
            
            for (m in 1:length(miss_vec)) {
                missr <- miss_vec[m]
                
                print(paste("Currently at n, or, missingness:", n_vec[j], or_vec[i], miss_vec[m]))
                
                # calculate row-major order index calculation (1-based)
                # from Wikipedia: https://en.wikipedia.org/wiki/Row-_and_column-major_order
                row_index <- ((j - 1) * length(or_vec) * length(miss_vec)) +
                    +     ((i - 1) * length(miss_vec)) + m
                
                # load result dfs from successful replications
                fname <- paste('df_results_succ',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                df_results_succ <- get(load(fname) )  # get object of loaded file
                
                ## create df with summarized results
                df_summ_results <- get_summary_stats(df_results_succ, n, or, missr, nreps)
                
                fname <- paste('df_summ_results',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                vname <- paste('df_summ_results',n,'_',or,'_',missr,impmethod, sep='')
                assign(vname, df_summ_results)
                
                
                # save result dfs as files
                print(paste("Saving file:", vname))
                save(df_summ_results, file = paste0(filepath_results, fname), compress = F)
                
                # save result row at correct index of result dataframe
                df_summ_all_cond[row_index, result_vars] <- df_summ_results

            }
        }
    }
    print("Saving final result file:")
    print(df_summ_all_cond)
    save(df_summ_all_cond, file = paste0(filepath_results, "df_summ_all_cond.RData"), compress = F)
}






