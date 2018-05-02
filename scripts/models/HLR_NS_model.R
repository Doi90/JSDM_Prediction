#########################################################
#########################################################    
###                                                   ###
###                  RUN HMSC JSDM                    ###    
###                                                   ###
### This script runs a HMSC JSDM set up for           ###
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

y <- as.matrix(y)                                       # Required format

## Site data

command <- sprintf("read.csv('data/%1$s/X_%1$s_fold%2$s_train_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

X <- eval(parse(text = command))                        # Evaluate command to read in data

X <- X[ , -1]                                           # Remove rownames

X <- as.matrix(X)                                       # Required format

## Set up latent factors

LF <- as.data.frame(rep("a",                            # Need to define a sampling level
                        nrow(X)))                       #  Here all same level

### Run HMSC JSDM ----

## Set MCMC parameters

n_iter <- 50000

n_burn <- 10000

n_thin <- 50

n_samples <- (n_iter - n_burn) / n_thin

## Format Data

form_data <- as.HMSCdata(Y = y,                 # Pres/Abs data
                         X = X,                 # Environmental data
                         Random = LF,           # Latent factor sampling levels
                         scaleX = FALSE,        # Don't standardise (done a priori)
                         interceptX = TRUE)     # Add an intercept column

form_prior <- as.HMSCprior(form_data,           # Feed in data
                           family = "probit")   # Generate priors for probit model

form_param <- as.HMSCparam(form_data,           # Feed in data
                           form_prior)          # Feed in priors

set.seed(28041948)                              # Creator's birthday

JSDM <- hmsc(data = form_data,                  # Formatted data
             priors = form_prior,               # Formatted priors
             param = form_param,                # Formatted parameters
             family = "probit",                 # Probit link for P/A data
             niter = n_iter,                    # MCMC iterations
             nburn = n_burn,                    # MCMC burn in
             thin = n_thin,                     # MCMC thinning
             verbose = FALSE)                   # Don't print progress updates to screen

### Extract posteriors ----

## Beta

Beta_extract <- JSDM$results$estimation$paramX  # Extract Beta posterior from model (wrong shape)

Beta_posterior <- array(NA,                     # Create empty array of required shape
                        dim = c(ncol(X),
                                ncol(y),
                                n_samples))

for(i in seq_len(n_samples)){                        # Extract samples and fill Beta_posterior
                                                     #  i = sample
  for(j in seq_len(length(Beta_extract))){           # Fill correct shape with posterior values                 
                                                     #  j = species
    Beta_posterior[ , j, i] <- Beta_extract[j, , i]  # Fill Beta_posterior
    
  }
}

## R

R_extract <- JSDM$results$estimation$paramLatent # Extract factor loadings posterior
                                                 #  from model (wrong shape)

R_posterior <- array(NA,                         # Create empty array of required shape
                     dim = c(ncol(y),
                             ncol(y),
                             nrow(n_samples)))

for(i in seq_len(dim(R_posterior)[3])){    # Fill correct shape with posterior values
  
  tmp <- R_extract[[i]]                    # Extract sample's factor loadings matrix
  
  tmp <- tmp %*% t(tmp)                    # Factor loadings => covariance matrix
  
  tmp <- tmp + diag(ncol(y))               # Need to add a diagonal matrix
  
  tmp <- cov2cor(tmp)                      # Convert covariance to correlation
  
  R_posterior[ , , i] <- tmp               # Save each full correlation matrix to an 
                                           #  array slice
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