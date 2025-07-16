
#--------------------------------------------------------------------------------
## Title:                       Functions for a Binary Simulation Study
## Author:                      Ermioni Athanasiadi
## Last modification:           10.07.2025
## binary outcome and simple one-factor experimental design
## manipulated conditions:
## manipulated conditions:
##      -- sample size n
##      -- effect size (coefficients: or, interc)
#--------------------------------------------------------------------------------

# this does not run from the console  # set working directory to where file is located
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#--------------------------------------------------------------------------------

check_even <- function(n) {
    if (any(n %% 2 != 0)) {
        stop("Please select even number as sample size")
    } 
}


#--------------------------------------------------------------------------------
    
library(dplyr)
library(ggplot2)
library(mice)
library(future.apply)
    
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

#--------------------------------------------------------------------------------
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

#--------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------
# Function to impute missing values given a dataframe and a method
impute_method <- function(dat, method) {
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

#-------------------------------------------------------------------------------

# performance_results:

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



#-------------------------------------------------------------------------------
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
                
                # generate imputated datasets for all replications
                time.parallel <- system.time(
                    list_sims <- future_lapply(
                        list_sims, function(i) impute_method(
                            i, 
                            method = impmethod),
                        future.seed = TRUE))
                
                # print execution time for imputation
                print(time.parallel)
                
                ## set filepath for result folder
                filepath = "./../../binary_simulation_study/results/"
                
                # save convergence information as df to file ################################
                ######
                # future idea:
                # only variable to be retrieved in this step, all others later
                # --> rename for clarity ("convergence"), save as vectors?
                ######
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
                
            }
        }
    }
}


#-------------------------------------------------------------------------------
# function to estimate effects for single dataframe given a list of imputed datasets
estimate_effects <- function(imp) {
    
    # initiate empty results df
    model_coefs <- data.frame("estimate" = NA, "std.error" = NA, "statistic" = NA, "df" = NA, "p.value" = NA)
    
        # only estimate pooled results if replication is successful
        if (!is.null(imp)) {
            imp_models <- with(data = imp, glm(y_i ~ group, family = "binomial"))
            model_coefs <- summary(mice::pool(imp_models))[-1, -1]
        }

    return(model_coefs)
}

#-------------------------------------------------------------------------------
# Function to calculate performance-related outcomes, given a dataframe
add_perf_results <- function(df)  { #what other parameters? 
    # given: a dataframe containing all rows with successful replications
    
    # boolean: significant p-value y=1 or n=0
    df$signif <- ifelse(df$p.value < .05, 1, 0)
    
    # Bias
    # abs_bias for or
    # rel_bias for or
    
    # abs_bias for se
    # rel_bias for se
    
    # percent zeros
    
    # runtimes ##################
    
    
    return(df)
}

#-------------------------------------------------------------------------------
# Function to analyse data and export dfs with results for all replications and conditions 
analyse_data <- function(n_vec, or_vec, interc, nreps, miss_vec, impmethod, seed) {
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
                
                # estimate effects for all replications (NA rows for non-converged replications)
                time.parallel <- system.time(
                    estimated_effects <- future_lapply(list_sims, function(i) estimate_effects(i), future.seed = TRUE)
                )
                
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
                # set file names
                fname <- paste('df_results_all',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                vname <- paste('df_results_all',n,'_',or,'_',missr,impmethod, sep='')
                assign(vname, df_results_all)
                # save complete result dfs as files
                print(paste("Saving file:", vname))
                save(df_results_all, file = paste0(filepath_results, fname), compress = F)
                
                # now look at subset of successful replications (convergence)
                # extract relevant rows 
                df_results_succ <- df_results_all %>%
                    filter(convergence == 1)
                
                # add further variables for performance-related outcomes #
                df_results_succ <- add_perf_results(df_results_succ)
                
                # set file names
                fname <- paste('df_results_succ',n,'_',or,'_',missr,impmethod,'.RData', sep='')
                vname <- paste('df_results_succ',n,'_',or,'_',missr,impmethod, sep='')
                assign(vname, df_results_succ)
                # save successful result dfs as files
                print(paste("Saving file:", vname))
                save(df_results_succ, file = paste0(filepath_results, fname), compress = F)

                
            }
        }
    }
}

#-------------------------------------------------------------------------------
# Function to summarize the results of successful replications of a simulation study, given a df
get_summary_stats <- function(df, n, or, missr) {
    # input: a dataframe with results for all successful replications, other pars
    
    mean_or_hat   <- mean(exp(df[, 1]))

    mean_prob_hat <- mean(exp(df[, 1])) / (1 + mean(exp(df[, 1])))
    power         <- sum(df$signif) / nrow(df)
    true_or       <- or
    true_prob     <- or / (1 + or)   # check this!!
    # print(paste("printing mean estimated or, mean_prob_hat:", mean_or_hat, mean_prob_hat))
    
    # absolute bias
    # relative bias
    
    # mean_se_hat
    # true_se
    
    # perc_zeros
    
    df_summ <- data.frame(mean_or_hat, mean_prob_hat, power, true_or, true_prob)

    return(df_summ)
    
}



#-------------------------------------------------------------------------------
# Function to summarize and export results as dfs for all simulation conditions
summarize_data <- function(n_vec, or_vec, interc, nreps, miss_vec, impmethod_vec, seed) {

    # extend this later to run for all impmethods ? ###########################
    impmethod <- impmethod_vec[1]
    ###########################################################################
    
    # specify filepath for retrieving and saving results
    filepath_results = "./../results/"
    fname <- paste('df_results_succ',n_vec[1],'_',or_vec[1],'_',miss_vec[1],impmethod,'.RData', sep='')
    df_results_succ <- get(load(fname) )  # get object of loaded file
    
    # initialize empty result dataframe for all conditions (except impmethod, see above)
    df_summ_all_cond <- expand.grid(impmethod = impmethod, 
                                    miss = miss_vec, 
                                    or = or_vec, 
                                    n = n_vec)[, c("n", "or", "miss", "impmethod")]
    # grab all result variables 
    result_vars <- names(get_summary_stats(df = df_results_succ, n = n_vec[1], or_vec[1], miss_vec[1]))
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
                df_summ_results <- get_summary_stats(df_results_succ, n, or, missr)
                
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





