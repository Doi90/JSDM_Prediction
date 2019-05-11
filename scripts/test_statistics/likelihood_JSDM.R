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

## Pres/Abs data

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_test.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y <- eval(parse(text = command))                        # Evaluate command to read in data

y <- y[ , -1]                                           # Remove rownames

## Site data

command <- sprintf("read.csv('data/%1$s/X_%1$s_fold%2$s_test.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

X <- eval(parse(text = command))                        # Evaluate command to read in data

X <- X[ , -1]                                           # Remove rownames

X <- cbind(1, X)

## Marginal probability prediction

filename <- sprintf("outputs/predictions/%s_%s_fold%s_marginal_prob.rds",
                    model_id,
                    dataset_id,
                    fold_id)

prediction <- readRDS(filename)

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

## Independent likelihood

independent_LL <- tryCatch(expr = independent_log_likelihood(y = y,
                                                             pred = prediction,
                                                             n_species = n_species,
                                                             n_sites = n_sites,
                                                             n_iter = n_iter),
                           error = function(e){ return(NA) })
## Save to file

filename <- sprintf("outputs/likelihood/%s_%s_fold%s_independent_likelihood.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(independent_LL,
        filename)

## Joint likelihood

joint_LL <- tryCatch(expr = joint_log_likelihood(Beta = Beta_posterior,
                                                 X = X,
                                                 y = y,
                                                 R = R_posterior,
                                                 n_species = n_species,
                                                 n_sites = n_sites,
                                                 n_iter = n_iter),
                     error = function(e){ return(NA) })
## Save to file

filename <- sprintf("outputs/likelihood/%s_%s_fold%s_joint_likelihood.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(joint_LL,
        filename)