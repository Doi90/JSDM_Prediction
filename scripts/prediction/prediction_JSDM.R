#############################################################
#############################################################
###                                                       ### 
###              DEFINE PREDICTION CODE                   ###
###                                                       ###
### This script defines the different types of prediction ###
### we use in our JSDM prediction comparison.             ###
###                                                       ###
#############################################################
#############################################################

#######################
### Load Posteriors ###
#######################

Beta_filename <- sprintf("posteriors/%s_beta_%s_fold_%s.rds",
                         model_id,
                         dataset_id,
                         fold_id)

Beta_posterior <- readRDS(Beta_filename)

R_filename <- sprintf("posteriors/%s_R_%s_fold_%s.rds",
                      model_id,
                      dataset_id,
                      fold_id)

R_posterior <- readRDS(R_filename)

#################
### Load Data ###
#################

### Presence/Absence ----

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_train_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y_train <- eval(parse(text = command))                  # Evaluate command to read in data

y_train <- y_train[ , -1]                               # Remove rownames

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_test_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y_test <- eval(parse(text = command))                   # Evaluate command to read in data

y_test <- y_test[ , -1]                                 # Remove rownames

### Environmental variables ----

# MUST INCLUDE INTERCEPTS!

command <- sprintf("read.csv('data/%1$s/X_%1$s_fold%2$s_train_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

X_train <- eval(parse(text = command))                  # Evaluate command to read in data

X_train <- X_train[ , -1]                               # Remove rownames

Intercept <- rep(1, nrow(X_train))                      # Create an intercept column

X_train <- cbind(Intercept, X_train)                    # Add intercept column to front of dataset

command <- sprintf("read.csv('data/%1$s/X_%1$s_fold%2$s_test_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

X_test <- eval(parse(text = command))                   # Evaluate command to read in data

X_test <- X_test[ , -1]                                 # Remove rownames

Intercept <- rep(1, nrow(X_test))                       # Create an intercept column

X_test <- cbind(Intercept, X_test)                      # Add intercept column to front of dataset
  
#----

#######################
### Load Posteriors ###
#######################

## Beta
# Array of regression coefficients.
# Row = Measured variable. K rows (INCLUDES INTERCEPT!)
# Column = Species. J columns
# Slice = Iterations. n_iter slices

Beta_RDS <- sprintf("readRDS('posteriors/%s_beta_%s_fold_%s.rds')",
                    model_id,
                    dataset_id,
                    fold_id)

Beta_posterior <- eval(parse(text = Beta_RDS))


## R
# Array of correlation coefficients.
# Row = Species. J rows
# Column = Species. J columns
# Slice = Iterations. n_iter slices

R_RDS <- sprintf("posteriors/%s_R_%s_fold_%s.rds",
                 model_id,
                 dataset_id,
                 fold_id)

R_posterior <- eval(parse(text=R_RDS))

########################
### Define Constants ###
########################

n_species <- ncol(y_train)
  
n_sites <- nrow(y_train)
  
n_covar <- ncol(X_train) # includes intercept

n_iter <-  dim(Beta_posterior)[3]  # Number of MCMC iterations in posterior chains. Post thinning
  


  


##################
### PREDICTION ###
##################

# DO NOT ALTER CODE IN THIS SECTION ONCE FINALISED
# MUST BE KEPT CONSTANT ACROSS ALL JSDMS FOR THE
# COMPARISON TO BE FAIR



### Run prediction code ----

marg_pred <- predict.marginal()

cond_LOI_pred <- predict.conditional.LOI()

cond_LOO_pred <- predict.conditional.LOO()

joint_pred <- predict.joint()

#----

### Evaluate test statistics ----


#----

##############
### OUTPUT ###
##############

prediction_list <- list(marg_pred = marg_pred,
                        marg_ts = marg_ts,
                        cond_LOI_pred = cond_LOI_pred,
                        cond_LOI_ts = cond_LOI_ts,
                        cond_LOO_pred = cond_LOO_pred,
                        cond_LOO_ts = cond_LOO_ts,
                        joint_pred = joint_pred,
                        joint_ts = joint_ts,
                        fold_id = fold_id,
                        dataset_id = dataset_id,
                        session_info = sessionInfo(),
                        date_time = Sys.time())
