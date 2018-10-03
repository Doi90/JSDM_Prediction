######################################################
######################################################
###                                                ###
###         JSDM LOG LIKELIHOOD WORKFLOW           ###
###                                                ###
###    This script calculates the log-likelihood   ###
### statistics used in this analysis for the JSDMs ###
###                                                ###
######################################################
######################################################

#######################
### Load Posteriors ###
#######################

## Beta
# Array of regression coefficients.
# Row = Measured variable. K rows (INCLUDES INTERCEPT!)
# Column = Species. J columns
# Slice = Iterations. n_iter slices

Beta_filename <- sprintf("outputs/posteriors/%s_beta_%s_fold_%s.rds",
                         model_id,
                         dataset_id,
                         fold_id)

Beta_posterior <- readRDS(Beta_filename)

## R
# Array of correlation coefficients.
# Row = Species. J rows
# Column = Species. J columns
# Slice = Iterations. n_iter slices

R_filename <- sprintf("outputs/posteriors/%s_R_%s_fold_%s.rds",
                      model_id,
                      dataset_id,
                      fold_id)

R_posterior <- readRDS(R_filename)

#################
### Load Data ###
#################

## Likelihood is calculated using the data used to fit the model
## so we need to load in the TRAINING data

## Pres/Abs data

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_train_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y <- eval(parse(text = command))                        # Evaluate command to read in data

y <- y[ , -1]                                           # Remove rownames

## Site data

command <- sprintf("read.csv('data/%1$s/X_%1$s_fold%2$s_train_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

X <- eval(parse(text = command))                        # Evaluate command to read in data

X <- X[ , -1]                                           # Remove rownames

########################
### Define Constants ###
########################

n_species <- ncol(y)

n_sites <- nrow(y)

n_covar <- ncol(X) # includes intercept

n_iter <-  dim(Beta_posterior)[3]  # Number of MCMC iterations in posterior chains. Post thinning

################################
### Calculate Log-Likelihood ###
################################

log_likelihood_estimate <- log_likelihood(Beta = Beta_posterior,
                                          X = X,
                                          y = y,
                                          R = R_posterior,
                                          n_species = n_species,
                                          n_sites = n_sites,
                                          n_iter = n_iter)

## Save to file

filename <- sprintf("outputs/likelihood/%s_%s_fold%s_likelihood.rds")

saveRDS(log_likelihood_estimate,
        filename)