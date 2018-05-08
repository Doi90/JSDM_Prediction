#############################################################
#############################################################
###                                                       ### 
###              PERFORM SSDM PREDICTION                  ###
###                                                       ###
### This script performs the SSDM model predictions.      ###
###                                                       ###
#############################################################
#############################################################

#######################
### Load Model List ###
#######################

filename <- sprintf("posteriors/%s_%s_fold_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

model_list <- readRDS(filename)

#################
### Load Data ###
#################

### Presence/Absence ----

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_test_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y_test <- eval(parse(text = command))                   # Evaluate command to read in data

y_test <- y_test[ , -1]                                 # Remove rownames

### Environmental variables ----

# MUST INCLUDE INTERCEPTS!

command <- sprintf("read.csv('data/%1$s/X_%1$s_fold%2$s_test_spatial.csv')", 
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

##################
### Prediction ###
##################

## Empty array

pred_array <- array(NA,
                    dim = c(n_sites,
                            n_species,
                            1))

## Predict for each species

for(j in seq_len(n_species)){
  
  new_data <- cbind(y_test[ , j],                       # Create single species dataframe
                    X_test)
  
  command <- sprintf("pred_array[ , %1$s, ] <-
                         predict.glm(object = model_list[[%1$s]],
                                     newdata = new_data,
                                     type = 'response')",
                     j,
                     colnames(y_test)[j])
  
  eval(parse(text = command))                           # Run command
  
}

###################
### Save Output ###
###################

filename <- sprintf("outputs/predictions/%s_%s_fold%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)
saveRDS(pred_array,
        filename)
