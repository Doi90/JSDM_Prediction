#########################################################
#########################################################    
###                                                   ###
###               RUN BayesComm JSDM                  ###    
###                                                   ###
### This script runs a BayesComm JSDM set up for      ###
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

### Run BayesComm JSDM ----

set.seed(28041948)                  # Creator's Birthday

JSDM <- BC(Y = y,                   # Pres/Abs data
           X = X,                   # Predictors
           model = "full",          # Coef+Corr model
           its = 11000,             # Number of MCMC samples
           burn = 1000,             # Size of burn-in
           thin = 10)               # Amount of thinning

### Extract posteriors ----

## Beta

Beta_extract <- JSDM$trace$B                    # Extract Beta posterior from model (wrong shape)

n_samples <- (JSDM$call$its - JSDM$call$start + 1) / JSDM$call$thin  # determine number of samples in posterior

Beta_posterior <- array(NA,                     # Create empty array of required shape
                        dim = c(ncol(X),
                                ncol(y),
                                n_samples))

for(j in seq_len(length(Beta_extract))){        # Fill correct shape with posterior values
                                                # j = species
  tmp <- Beta_extract[[j]]                      # Extract single species from list
  
  for(i in seq_len(n_samples)){                 # Extract samples and fill Beta_posterior
                                                # i = sample
    Beta_posterior[ , j, i] <- tmp[i, ]         # Fill Beta_posterior
    
  }
}

## R

R_extract <- JSDM$trace$R                       # Extract R posterior from model (wrong shape)
    
R_posterior <- array(NA,                        # Create empty array of required shape
                     dim = c(ncol(y),
                             ncol(y),
                             nrow(R_extract)))

for(i in seq_len(dim(R_posterior)[3])){         # Fill correct shape with posterior values
  
  tmp <- matrix(0,                              # 0 not NA so we can use addition later
                nrow = ncol(y),
                ncol = ncol(y))
  
  tmp[upper.tri(tmp)] <- R_extract[i, ]         # Fill upper.tri of matrix
  
  tmp <- tmp + t(tmp)                           # Fill lower.tri of matrix
  
  diag(tmp) <- 1                                # Ones on diagonals
  
  R_posterior[ , , i] <- tmp                    # Save each full correlation matrix to an 
                                                # array slice
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