
## Script to generate data for a simulation study
## binary outcome and simple one-factor experimental design
## manipulated conditions:
##      -- sample size n
##      -- effect size (coefficients: or, interc)

# this does not run from the console
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


# helper Function to generate Data for single dataset
generate_df <- function(n, or, interc) {
    
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

# Function to create replications of generated datasets
generate_replications <- function(n_vec, or_vec, interc, nreps, seed) {
    set.seed(seed)
    for (j in 1:length(n_vec)) {
        n <- n_vec[j]
        
        for (i in 1:length(or_vec)) {
            or <- or_vec[i]
            time.parallel <- system.time(list_sims <- future_lapply(1:nreps, function(i) generate_df(n = n, or = or, interc = interc), future.seed = TRUE))
            fname <- paste('list_data',n,'_',or,'_',j,i,'.RData', sep='')
            save(list_sims, file = fname, compress = F)
            
        }
    }
}


# export execution time ###
# time.parallel
