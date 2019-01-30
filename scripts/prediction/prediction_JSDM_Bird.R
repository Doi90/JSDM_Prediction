#############################################################
#############################################################
###                                                       ### 
###              PERFORM BIRD PREDICTION                  ###
###                                                       ###
### This script performs the types of prediction defined  ###
### for JSDMs in this comparison. Performs each type of   ###
### prediction in series, & saves to file. This script    ###
### is ONLY for the bird dataset runs as it involves      ###
### running on a subset of posterior samples (10) only.   ###
###                                                       ###
#############################################################
#############################################################

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

Beta_posterior <- Beta_posterior[ , , start_sample:end_sample]

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

R_posterior <- R_posterior[ , , start_sample:end_sample]

#################
### Load Data ###
#################

### Presence/Absence ----

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_test.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y_test <- eval(parse(text = command))                   # Evaluate command to read in data

y_test <- y_test[ , -1]                                 # Remove rownames

### Environmental variables ----

# MUST INCLUDE INTERCEPTS!

command <- sprintf("read.csv('data/%1$s/X_%1$s_fold%2$s_test.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

X_test <- eval(parse(text = command))                   # Evaluate command to read in data

X_test <- X_test[ , -1]                                 # Remove rownames

Intercept <- rep(1, nrow(X_test))                       # Create an intercept column

X_test <- cbind(Intercept, X_test)                      # Add intercept column to front of dataset

#----

########################
### Define Constants ###
########################

n_species <- ncol(y_test)

n_sites <- nrow(y_test)

n_covar <- ncol(X_test) # includes intercept

n_iter <-  dim(Beta_posterior)[3]  # Number of MCMC iterations in posterior chains. Post thinning

##################
### Prediction ###
##################

## Marginal

marg_start <- Sys.time()

## Probabilities

marg_pred_prob <- predict.marginal.probability(Beta = Beta_posterior,
                                               X = X_test,
                                               n_species = n_species,
                                               n_sites = n_sites,
                                               n_iter = n_iter)

filename <- sprintf("outputs/predictions/%s_%s_fold%s_marginal_prob_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(marg_pred_prob,
        filename)

## Binary

marg_pred_bin <- predict.marginal.binary(Beta = Beta_posterior,
                                         X = X_test,
                                         n_species = n_species,
                                         n_sites = n_sites,
                                         n_iter = n_iter)

filename <- sprintf("outputs/predictions/%s_%s_fold%s_marginal_bin_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(marg_pred_bin,
        filename)

message(sprintf("Marginal prediction duration: %s hours",
                round(difftime(Sys.time(),
                               marg_start,
                               units = "hours")[[1]],
                      digits = 5)))

rm(marg_start,
   marg_pred_prob,
   marg_pred_bin)

# ## Conditional - Leave One Out
# 
# cond_LOO_start <- Sys.time()
# 
# cond_LOO_pred <- predict.conditional.LOO(Beta = Beta_posterior,
#                                          X = X_test,
#                                          y = y_test,
#                                          R = R_posterior,
#                                          n_species = n_species,
#                                          n_sites = n_sites,
#                                          n_iter = n_iter)
# 
# filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOO_%s_%s.rds",
#                     model_id,
#                     dataset_id,
#                     fold_id,
#                     start_sample,
#                     end_sample)
# 
# saveRDS(cond_LOO_pred,
#         filename)
# 
# message(sprintf("Conditional LOO prediction duration: %s hours",
#                 round(difftime(Sys.time(),
#                                cond_LOO_start,
#                                units = "hours")[[1]],
#                       digits = 5)))
# 
# rm(cond_LOO_start,
#    cond_LOO_pred)

## Conditional - Leave One In

cond_LOI_start <- Sys.time()

cond_LOI_pred <- predict.conditional.LOI(Beta = Beta_posterior,
                                         X = X_test,
                                         y = y_test,
                                         R = R_posterior,
                                         n_species = n_species,
                                         n_sites = n_sites,
                                         n_iter = n_iter,
                                         dataset_id = dataset_id)

filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOI_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(cond_LOI_pred,
        filename)

message(sprintf("Conditional LOI prediction duration: %s hours",
                round(difftime(Sys.time(),
                               cond_LOI_start,
                               units = "hours")[[1]],
                      digits = 5)))

rm(cond_LOI_start,
   cond_LOI_pred)

## Joint

joint_start <- Sys.time()

joint_pred <- predict.joint(Beta = Beta_posterior,
                            X = X_test,
                            y = y_test,
                            R = R_posterior,
                            n_species = n_species,
                            n_sites = n_sites,
                            n_iter = n_iter)

filename <- sprintf("outputs/predictions/%s_%s_fold%s_joint_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(joint_pred,
        filename)

message(sprintf("Joint prediction duration: %s hours",
                round(difftime(Sys.time(),
                               joint_start,
                               units = "hours")[[1]],
                      digits = 5)))

rm(joint_start,
   joint_pred)
