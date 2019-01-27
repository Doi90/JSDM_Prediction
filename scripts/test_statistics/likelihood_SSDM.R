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

SSDM_predictions <- readRDS(filename)

## SESAM

filename <- sprintf("outputs/predictions/SESAM_%s_fold%s.rds",
                    dataset_id,
                    fold_id)

SESAM_predictions <- readRDS(filename)

################################
### Calculate Log-Likelihood ###
################################ 

# SSDM 

## Vector to store likelihoods

SSDM_log_lik_vector <- vector(length = dim(SSDM_predictions)[1])

## Calculate likelihoods

for(i in seq_len(dim(SSDM_predictions)[1])){
  
  ## Get predicted species probabilities for site
  
  pred_probs <- SSDM_predictions[i, , 1]
  
  ## Get observed species state for site
  
  obs_spp <- y[i, ]
  
  ## Calculate log likelihood for site
  
  log_lik <- sum(dbinom(x = unlist(obs_spp),
                        size = 1,
                        prob = pred_probs,
                        log = TRUE))
  
  ## Add value to vector
  
  SSDM_log_lik_vector[i] <- log_lik
  
}

## Obtain single likelihood value for whole model

SSDM_log_likelihood_estimate <- sum(SSDM_log_lik_vector,
                                    na.rm = TRUE)

## Save to file

filename <- sprintf("outputs/likelihood/SSDM_%s_fold%s_likelihood.rds",
                    dataset_id,
                    fold_id)

saveRDS(SSDM_log_likelihood_estimate,
        filename)

# SESAM

## Vector to store likelihoods

SESAM_log_lik_vector <- vector(length = dim(SESAM_predictions)[1])

## Calculate likelihoods

for(i in seq_len(dim(SESAM_predictions)[1])){
  
  ## Get predicted species probabilities for site
  
  pred_probs <- SESAM_predictions[i, , 1]
  
  ## Get observed species state for site
  
  obs_spp <- y[i, ]
  
  ## Calculate log likelihood for site
  
  log_lik <- sum(dbinom(x = unlist(obs_spp),
                        size = 1,
                        prob = unlist(pred_probs),
                        log = TRUE))
  
  ## Add value to vector
  
  SESAM_log_lik_vector[i] <- log_lik
  
}

## Obtain single likelihood value for whole model

SESAM_log_likelihood_estimate <- sum(SESAM_log_lik_vector,
                                     na.rm = TRUE)

## Save to file

filename <- sprintf("outputs/likelihood/SESAM_%s_fold%s_likelihood.rds",
                    dataset_id,
                    fold_id)

saveRDS(SESAM_log_likelihood_estimate,
        filename)
