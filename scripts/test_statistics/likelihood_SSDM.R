######################################################
######################################################
###                                                ###
###      SSDM/SESAM LOG LIKELIHOOD WORKFLOW        ###
###                                                ###
###    This script calculates the log-likelihood   ###
### statistics used in this analysis for SSDMs and ###
### SESAM models.                                  ###
###                                                ###
######################################################
######################################################

#################
### Load Data ###
#################

## Likelihood is calculated using the data used to fit the model
## so we need to load in the TRAINING data

## Pres/Abs data

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_train.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y <- eval(parse(text = command))                        # Evaluate command to read in data

y <- y[ , -1]                                           # Remove rownames

########################
### Load Predictions ###
########################

## SSDM

filename <- sprintf("outputs/predictions/%s_%s_fold%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

readRDS(filename)

## SESAM

filename <- sprintf("outputs/predictions/SESAM_%s_fold%s.rds",
                    dataset_id,
                    fold_id)

readRDS(filename)

################################
### Calculate Log-Likelihood ###
################################ 

## Vector to store likelihoods

log_lik_vector <- vector(length = dim(pred_array)[1])

## Calculate likelihoods

for(i in seq_len(dim(pred_array)[1])){
  
  ## Get predicted species probabilities for site
  
  pred_probs <- pred_array[i, , 1]
  
  ## Get observed species state for site
  
  obs_spp <- y[i, ]
  
  ## Calculate log likelihood for site
  
  log_lik <- sum(dbinom(x = obs_spp,
                        size = 1,
                        prob = pred_probs,
                        log = TRUE))
  
  ## Add value to vector
  
  log_lik_vector[i] <- log_lik
  
}

## Obtain single likelihood value for whole model

log_likelihood_estimate <- sum(log_lik_vector,
                               na.rm = TRUE)

## Save to file

filename <- sprintf("outputs/likelihood/%s_%s_fold%s_likelihood.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(log_likelihood_estimate,
        filename)
