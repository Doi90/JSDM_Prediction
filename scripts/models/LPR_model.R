#########################################################
#########################################################    
###                                                   ###
###                  RUN boral JSDM                   ###    
###                                                   ###
### This script runs a boral JSDM set up for          ###
### cross-validation on Spartan. Also defines outputs ###
### in correct format for prediction code.            ###  
###                                                   ###
#########################################################
#########################################################

### Load Data ----

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

### Run boral JSDM ----

JSDM <- boral(y,                                           # Pres/Abs data
              X = X,                                       # Predictors. DO NOT INCLUDE INTERCEPT COLUMN
              family = "binomial",                         # Makes model use PA data, probit link
              num.lv = 2,                                  # Set number of latent factors
              save.model = TRUE,                           # Saves JAGS model as a txt file, allows coda package to analyse MCMC
              mcmc.control = list(n.burnin = 10000,        # Burn in
                                  n.iteration = 60000,     # Total number of samples
                                  n.thin = 50,             # Amount of thinning
                                  seed = 28041948),        # Seed
              model.name = sprintf("boral_JAGS_%s_%s.txt", # Name of saved txt file.
                                   dataset_id,
                                   fold_id))                         

### Extract posteriors ----

## Beta

Beta_extract <- JSDM$jags.model$BUGSoutput$sims.list    # Extract Beta posterior from model (wrong shape)

n_samples <- (JSDM$mcmc.control$n.iteration -
                JSDM$mcmc.control$n.burnin) / JSDM$call$thin  # Determine number of samples in posterior

Beta_posterior <- array(NA,                     # Create empty array of required shape
                        dim = c(ncol(X),
                                ncol(y),
                                n_samples))

for(j in seq_len(ncol(Beta_extract$X.coefs))){  # Fill correct shape with posterior values
                                                #  j = species
  for(i in seq_len(n_samples)){                 # Extract samples and fill Beta_posterior
                                                #  i = sample
    tmp <- c(Beta_extract$lv.coefs[i, j, 1],    # Intercept
             Beta_extract$X.coefs[i, j, ])      # Non-intercept regression coefficients
    
    Beta_posterior[ , j, i] <- tmp              # Fill Beta_posterior
    
  }
}

## R

R_extract <- JSDM$jags.model$BUGSoutput$sims.list$lv.coefs  # Extract factor loadings posterior from model
R_extract <- R_extract[ , , -1]                             # Drop the Intercept values

R_posterior <- array(NA,                        # Create empty array of required shape
                     dim = c(ncol(y),
                             ncol(y),
                             n_samples))

for(i in seq_len(dim(R_posterior)[3])){         # Fill correct shape with posterior values
                                                #  i=samples
  lambda <- R_extract[i, , ]                    # Extract factor loadings matrix
                                                #  lambda = J rows, H col
  tmp <- lambda %*% t(lambda) + diag(ncol(y))   # Convert factor loadings to covariance matrix
  
  tmp <- cov2cor(tmp)                           # Convert covariance matrix to correlation matrix
  
  R_posterior[ , , i] <- tmp                    # Save each full correlation matrix to an array slice

}

### Save Posteriors ----

filename <- sprintf("posteriors/%s_beta_%s_fold_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(Beta_posterior,
        filename)

filename <- sprintf("posteriors/%s_R_%s_fold_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(R_posterior,
        filename)