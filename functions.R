
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
    
    ## add check for updates:
    # remotes::update_packages('rblimp')
    
    dat$y_i <- as.integer(dat$y_i)
    dat$y_i <- ifelse(dat$y_i == 1, 0, 1)
    
    # truemodel <- glm(y_i ~ group, data = dat, family = binomial())
    # summary(truemodel)
    
    # tryCatch to handle any errors and return NULL if imputation fails
    imp_blimp <- tryCatch({
        suppressWarnings({
            rblimp(
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
            as.mitml(imp_blimp)
        })
    }, error = function(e) {
        message("XXX Error in imputation: ", e$message)  # print the error message
        return(NULL)  # return NULL to indicate failure
    }) 
    
    # if (is.null(imp_blimp)) {
    #     return(NULL)  # return integer 0 explicitly
    # }

    return(imp_blimp)
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

                # if valid impmethod, generate imputed datasets for all replications
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
                    
                    # save execution time to list for later access
                    list_sims <- lapply(list_sims, function(x) {
                        if (!is.null(x)) {
                            attr(x, "timing")  <- time.parallel 
                            attr(x, "timing2") <- time.parallel[3]
                        }
                        x
                    })
                    
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
                
                # now add timing to df
                # move this into performance_results_function?
                df_perf_results$timing <- sapply(list_sims, function(x) attr(x, "timing")[3])
                
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
    # correction ?
    # df$convergence <- ifelse(sapply(list_sims, is.null), 0, 1)
    
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
    model <- glm(y_i ~ group, family = "binomial", data = df)
    model_summ <- summary(model)$coefficients[-1, ]
    
    # extract relevant parameters for output
    model_coefs$conf.low <- confint(model)[2]
    model_coefs$conf.high <- confint(model)[4]
    model_coefs$estimate <- model_summ[1]
    model_coefs$statistic <- NA # keep column because present in results from imp objects
    model_coefs[, c(1, 2, 5)] <- c(model_summ[1], model_summ[2], model_summ[4])
    
    return(model_coefs)
}


## UPDATED
estimate_effects <- function(df) {
    # initiate empty results df
    model_coefs <- data.frame(
        estimate   = NA,
        std.error  = NA,
        statistic  = NA,
        df         = NA,
        p.value    = NA,
        conf.low   = NA,
        conf.high  = NA
    )
    
    # quick checks before fitting
    if (nrow(df) == 0) {
        return(model_coefs)
    }
    
    # if group has fewer than 2 levels after LD → cannot fit model
    if (length(unique(df$group)) < 2) {
        return(model_coefs)
    }
    
    # fit logistic regression safely
    fit <- tryCatch(
        glm(y_i ~ group, family = "binomial", data = df),
        error = function(e) return(NULL),
        warning = function(w) invokeRestart("muffleWarning") # suppress warnings
    )
    
    # if fitting failed, return NA row
    if (is.null(fit)) {
        return(model_coefs)
    }
    
    # extract coefficients safely
    model_summ <- summary(fit)$coefficients
    # skip intercept
    if (nrow(model_summ) < 2) {
        return(model_coefs)
    }
    
    # confidence intervals (use tryCatch, because confint may fail too)
    ci <- tryCatch(confint(fit), error = function(e) matrix(NA, ncol = 2, nrow = 2))
    
    # fill in coefficients
    model_coefs$estimate  <- model_summ[2, 1]
    model_coefs$std.error <- model_summ[2, 2]
    model_coefs$p.value   <- model_summ[2, 4]
    model_coefs$statistic <- NA  # not available for glm by default
    model_coefs$conf.low  <- ci[2, 1]
    model_coefs$conf.high <- ci[2, 2]
    
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

analyse_complete_data <- function(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete) {
    
    for (j in 1:length(n_vec)) {
        n <- n_vec[j]
        
        for (i in 1:length(or_vec)) {
            or <- or_vec[i]
            

            print(getwd())
            print(paste("Currently at n, or:", n_vec[j], or_vec[i]))
            
            
            # if is.complete == TRUE
            if (impmethod == "none") {
                fname <- paste('list_data',n,'_',or,'.RData', sep='')
            }
            else {
                print("not yet coded")
                #fname <- paste('list_data_imp',n,'_',or,'_',missr,impmethod,'.RData', sep='')
            }
            
            
            print(fname)
            
            list_sims <- get(load(fname) )  # get object of loaded file

            # specify filepath for retrieving and saving results
            filepath_results = "./../../binary_simulation_study/results/"
            # retrieve df with performance results
            
            fname_perf <- paste('df_perf_results',n,'_',or,'_',missr,impmethod,'.RData', sep='')
            
            df_perf_results <- get(load(paste0(filepath_results, fname_perf)))
            
            
            #df_perf_results <- add_perf_results()
            
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
            
            ### add check here AND later: 
            ## only calculate the performance measures if at least 50% 
            # of the replications were successful
            
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


#________________________________________________________________________________
# Function to analyse data and export dfs with results for all replications and conditions 
analyse_data <- function(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed, is.complete) {
    
    for (j in 1:length(n_vec)) {
        n <- n_vec[j]
        
        for (i in 1:length(or_vec)) {
            or <- or_vec[i]
            
            # different missingness levels 
            for (m in 1:length(miss_vec)) {
                print(getwd())
                print(paste("Currently at n, or, missingness:", n_vec[j], or_vec[i], miss_vec[m]))
                
                # load imputed dataframes
                missr <- miss_vec[m]
                
                # if is.complete == TRUE
                # fname <- paste('list_data',n,'_',or,'.RData', sep='')
                
                if (impmethod == "none") {
                    # fname <- paste('list_data_imp',n,'_',or,'_',missr,'.RData', sep='')
                    fname <- paste('list_data',n,'_',or,'.RData', sep='')
                }
                
                else if (impmethod == "ld") {
                    fname <- paste('list_data_na',n,'_',or,'_',missr,'.RData', sep='')
                    # placeholder to retain same structure - no timing and convergence info for ld
                    df_perf_results <- data.frame(convergence = rep(NA, nreps), timing = rep(NA, nreps))
                }
                
                else {
                    fname <- paste('list_data_imp',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                    # retrieve df with performance results
                    fname_perf <- paste('df_perf_results',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                    df_perf_results <- get(load(paste0(filepath_results, fname_perf)))
                }
                

                print(fname)

                list_sims <- get(load(fname) )  # get object of loaded file
                # NOTE: depending on earlier condition check, this is either
                # a list of imputed datasets or complete or incomplete datasets
                # rewrite this later for easier debugging
 
                # specify filepath for retrieving and saving results
                filepath_results = "./../../binary_simulation_study/results/"


                
                # if no impmethod specified or data is complete: simple effect estimation
                # Note: this performs listwise deletion if data is incomplete
                if (impmethod == "none" | is.complete == TRUE) {
                    time.parallel <- system.time(
                    estimated_effects <- future_lapply(list_sims, function(i) estimate_effects(i), future.seed = TRUE)
                )
                }
                
                else if (impmethod == "ld") {
                    # listwise deletion
                    print(paste0("Listwise deletion: now estimating effects for ", fname))
                    time.parallel <- system.time(
                        estimated_effects <- future_lapply(
                            list_sims, # loads datasets with NAs
                            function(i) estimate_effects(i), 
                            future.seed = TRUE
                        )
                    )
                    df_perf_results$timing <- time.parallel[3]
                }
                
                # if impmethod specified: estimate effects with imp object
                else {
                    # estimate effects for all replications (NA rows for non-converged replications)
                    print(paste0("Imputation method: ", impmethod, " now estimating effects for ", fname))
                    time.parallel <- system.time(
                        estimated_effects <- future_lapply(
                            list_sims, # loads imputed datasets
                            function(i) estimate_effects_imp(i), 
                            future.seed = TRUE)
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
                
                ## ADD:
                # not only check imputation NULL objects, but also failed regression
                # analyses 

                
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
                    # print("correctly checking that we have a complete dataset here")
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
                    filter(convergence == 1 | is.na(convergence))  # add na condition
                                                                   # for ld case
                
                ### add check here AND later: 
                ## only calculate the performance measures if at least 50% 
                # of the replications were successful
                
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
get_summary_stats <- function(df, n, or, missr, nreps, flag_success) {
    # input: a dataframe with results for all successful replications, other pars
    #print(df)
    
    # if no
    if (flag_success == TRUE & nrow(df) != 0) {
        
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
    
    #runtime        <- NA
    ## runtimes:
    # no still given as col with repeated equal number: sum over all nreps
    if (!is.list(df$timing[1])) {
     runtime       <- df$timing[1] / nreps
    }
    
    #else if (is.null())
    
    else {
        print("timing is inside extra list element")
        runtime       <- df$timing[[1]] / nreps
    }

    
    #save(df$timing, file = "timing-test.RData")
    print(paste("current summary: n, or, missr", n, or, missr))
    #print(summary(nreps))
    print(df$timing[1])  # NA
    print(df$timing2)  # NULL
    
    runtime2 <- df$timing2[1]
    
    
    # perc_zeros
    df_summ <- data.frame(
        mean_or_hat, 
        mean_prob_hat, 
        true_or, 
        true_prob,
        mean_se_hat,
        abs_bias_or,
        rel_bias_or,
        power, 
        coverage,
        perc_convergence,
        runtime
    )
    
    }
    else {
        df_summ <- data.frame(
            mean_or_hat = NA, 
            mean_prob_hat = NA,
            true_or = NA,
            true_prob = NA,
            mean_se_hat = NA,
            power = NA,
            coverage = NA,
            perc_convergence = NA,
            runtime = NA
        )
    }
        


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
    result_vars <- names(get_summary_stats(df = df_results_succ, n = n_vec[1], or_vec[1], miss_vec[1], nreps, flag_success = TRUE))
    df_summ_all_cond[, result_vars] <- NA
    
    #print("TEST WORK UNTIL HERE?")
    
    # get results for each condition
    for (j in 1:length(n_vec)) {
        n <- n_vec[j]
        
        for (i in 1:length(or_vec)) {
            or <- or_vec[i]
            
            for (m in 1:length(miss_vec)) {
                missr <- miss_vec[m]
                
                print(paste("Currently at n, or, missingness:", n_vec[j], or_vec[i], miss_vec[m]))
                
                # calculate row-major order index (1-based)
                # from Wikipedia: https://en.wikipedia.org/wiki/Row-_and_column-major_order
                row_index <- ((j - 1) * length(or_vec) * length(miss_vec)) +
                    +     ((i - 1) * length(miss_vec)) + m
                
                # load result dfs from successful replications
                fname <- paste('df_results_succ',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                df_results_succ <- get(load(fname) )  # get object of loaded file
                
                ## add condition check e.g. ONLY estimate summary stats if at least 
                ## 30% (or 50%??) of replications were successful
                
                if (nrow(df_results_succ) >= nreps*0.5) {
                    print(paste0("Enough successful replications: ",
                                 nrow(df_results_succ), "out of", nreps, 
                                 "were successful. Obtaining summary statistics now"))
                    flag_success = TRUE
                }
                
                else {
                    print("flag works correctly:")
                    flag_success = FALSE
                }
                
                ## create df (single row) with summarized results
                ## if not enough successful runs, row with only NAs added
                df_summ_results <- get_summary_stats(df_results_succ, n, or, missr, nreps, flag_success)
                
                
                
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
    save(df_summ_all_cond, file = paste0(filepath_results, "df_summ_all_cond",impmethod, ".RData"), compress = F)
}


#--Part V: Visualization of Results -------------------------------------------------

#________________________________________________________________________________
# Function to visualize results, given a result dataframe and an imputation method

visualize_bias <- function(impmethod, show_title = TRUE, figwidth = 5, figheight = 3) {

    filepath_results = "./../results/"
    fname <- paste0(filepath_results, "df_summ_all_cond",impmethod, ".RData")
    
    df <- get(load(fname) )  # get object of loaded file
    
    if (show_title == TRUE) {
        plot_title = "Bias by Sample Size and Missingness Rate"
    }
    else
        plot_title = ""
    
    
    plt <- ggplot(df, aes(x = n, y = rel_bias_or, color = factor(true_or), group = true_or)) +
        geom_point() +  # size = 1.8
        geom_line() +  # linewith = 0.8
        facet_wrap(~ miss, labeller = label_both) +
        labs(
            title = plot_title,
            x = "n",
            y = "Bias (%)",
            color = "True OR"
        ) +
        theme_minimal()
    plt
    
    fpath = "./../results/figures/"
    type = "bias"
    fname <- paste0(fpath, type, "_", impmethod, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    
    ggsave(fname, plot = plt, width = figwidth, height = figheight)
    
}

visualize_power <- function(impmethod, show_title = TRUE, figwidth = 5, figheight = 3) {
    
    filepath_results = "./../results/"
    fname <- paste0(filepath_results, "df_summ_all_cond",impmethod, ".RData")
    df <- get(load(fname) )  # get object of loaded file
    
    if (show_title == TRUE) {
        plot_title = "Power by Sample Size and Missingness Rate"
    }
    else
        plot_title = ""
    
    plt <- ggplot(df, aes(x = n, y = power, color = factor(true_or), group = true_or)) +
        geom_point() +
        geom_line() +
        facet_wrap(~ miss, labeller = label_both) +
        labs(
            title = plot_title,
            x = "n",
            y = "Power",
            color = "True OR"
        ) +
        theme_minimal()
    plt
    fpath = "./../results/figures/"
    type = "power"
    fname <- paste0(fpath, type, "_", impmethod, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    
    ggsave(fname, plot = plt, width = figwidth, height = figheight)
    
    
}


visualize_coverage <- function(impmethod, show_title = TRUE, figwidth = 5, figheight = 3) {
    
    
    filepath_results = "./../results/"
    fname <- paste0(filepath_results, "df_summ_all_cond",impmethod, ".RData")
    df <- get(load(fname) )  # get object of loaded file
    
    if (show_title == TRUE) {
        plot_title = "Coverage by Sample Size and Missingness Rate"
    }
    else
        plot_title = ""
    
    plt <- ggplot(df, aes(x = n, y = coverage, color = factor(true_or), group = true_or)) +
        geom_point() +
        geom_line() +
        facet_wrap(~ miss, labeller = label_both) +
        labs(
            title = plot_title,
            x = "n",
            y = "Coverage",
            color = "True OR"
        ) +
        theme_minimal()
    
    # show figure
    plt
    
    fpath = "./../results/figures/"
    type = "coverage"
    fname <- paste0(fpath, type, "_", impmethod, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    
    ggsave(fname, plot = plt, width = figwidth, height = figheight)
    
    
}


visualize_convergence <- function(impmethod, show_title = TRUE, figwidth = 5, figheight = 3) {
    
    
    filepath_results = "./../results/"
    fname <- paste0(filepath_results, "df_summ_all_cond",impmethod, ".RData")
    df <- get(load(fname) )  # get object of loaded file
    
    if (show_title == TRUE) {
        plot_title = "Convergence by Sample Size and Missingness Rate"
    }
    else
        plot_title = ""
    
    plt <- ggplot(df, aes(x = n, y = perc_convergence, color = factor(true_or), group = true_or)) +
        geom_point() +
        geom_line() +
        facet_wrap(~ miss, labeller = label_both) +
        labs(
            title = plot_title,
            x = "n",
            y = "Convergence (%)",
            color = "True OR"
        ) +
        theme_minimal()
    
    plt
    
    fpath = "./../results/figures/"
    type = "convergence"
    fname <- paste0(fpath, type, "_", impmethod, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    
    ggsave(fname, plot = plt, width = figwidth, height = figheight)
    
    
}


visualize_runtime <- function(impmethod, show_title = TRUE, figwidth = 5, figheight = 3) {
    
    
    filepath_results = "./../results/"
    fname <- paste0(filepath_results, "df_summ_all_cond",impmethod, ".RData")
    df <- get(load(fname) )  # get object of loaded file
    
    if (show_title == TRUE) {
        plot_title = "Runtime by Sample Size and Missingness Rate"
    }
    else
        plot_title = ""
    
    plt <- ggplot(df, aes(x = n, y = runtime, color = factor(true_or), group = true_or)) +
        geom_point() +
        geom_line() +
        facet_wrap(~ miss, labeller = label_both) +
        labs(
            title = plot_title,
            x = "n",
            y = "Runtime",
            color = "True OR"
        ) +
        theme_minimal()
    
    plt
    
    fpath = "./../results/figures/"
    type = "runtime"
    fname <- paste0(fpath, type, "_", impmethod, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    
    ggsave(fname, plot = plt, width = figwidth, height = figheight)
    
}


visualize_runtime_by_impmethod <- function(impmethods = c("pmm", "logreg", "ld"),
                                           show_title = TRUE, 
                                           figwidth = 5, figheight = 3) {
    
    
    filepath_results = "./../results/"
    
    # dynamically load data for selected methods
    dfs <- lapply(impmethods, function(m) {
        fname <- paste0(filepath_results, "df_summ_all_cond", m, ".RData")
        df <- get(load(fname))
        df <- df %>% mutate(impmethod = m)
        return(df)
    })
    
    df_combined <- bind_rows(dfs)
    
    if (show_title == TRUE) {
        plot_title = "Runtime by Imputation Method, Sample Size and Missingness Rate"
    }
    else
        plot_title = ""
    
    
    plt <- ggplot(df_combined, aes(x = n, y = runtime, color = impmethod, group = impmethod)) +
        geom_point() +
        geom_line() +
        facet_wrap(~ miss, labeller = label_both) +
        labs(
            title = plot_title,
            x = "Sample Size (n)",
            y = "Runtime (s)",
            color = "Imputation Method"
        ) +
        theme_minimal()
    
    plt
    
    fpath = "./../results/figures/"
    type = "runtime_impmethod"
    fname <- paste0(fpath, type, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    
    
    ggsave(fname, plot = plt, width = figwidth, height = figheight)
    
    
    
}


get_all_visualizations <- function() {
    
    impmethods <- c("pmm", "logreg", "ld")
    
    lapply(impmethods, function(impmethod) {
    visualize_bias(impmethod, show_title = FALSE)
    visualize_power(impmethod, show_title = FALSE)
    visualize_coverage(impmethod, show_title = FALSE)
    visualize_convergence(impmethod, show_title = FALSE)
    visualize_runtime(impmethod, show_title = FALSE)
    }
    )
    visualize_runtime_by_impmethod(impmethods, show_title = FALSE)
}


