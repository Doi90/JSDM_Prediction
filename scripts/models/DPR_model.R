#########################################################
#########################################################    
###                                                   ###
###                  RUN gjam JSDM                    ###    
###                                                   ###
### This script runs a gjam JSDM set up for           ###
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

### Run gjam JSDM ----

form <- as.formula(paste("~ ", paste(colnames(X),       # Build formula from X column names
                                     collapse = "+")))

set.seed(28041948)                                      # Creator's birthday

JSDM <- gjam(formula = form,                              # Model formula
             xdata = X,                                   # Environment data
             ydata = y,                                   # Species data
             modelList = list(ng = 60000,                 # Total iterations
                              burnin = 10000,             # Burn-in
                              thin = 50,                  # Thinning
                              typeNames = 'PA',           # Data type = pres/abs
                              notStandard = colnames(X))) # Don't standardise

### Extract posteriors ----

## Beta

Beta_extract <- JSDM$chains$bgibbs                         # Extract Beta posterior from model (wrong shape)

Beta_extract <- Beta_extract[-(1:JSDM$modelList$burnin), ] # Remove burnin

Beta_extract <- Beta_extract[seq(1,                        # Implement thinning
                                 nrow(Beta_extract),
                                 JSDM$modelList$thin),] 

n_samples <- ((JSDM$modelList$ng - JSDM$modelList$burnin) / JSDM$modelList$thin)  # determine number of samples in posterior

Beta_posterior <- array(NA,                     # Create empty array of required shape
                        dim = c(ncol(X) + 1,    # Include intercept
                                ncol(y),
                                n_samples))

for(j in seq_len(ncol(y))){        # Fill correct shape with posterior values
                                                #  j = species
  array_id <- (j * (ncol(X) + 1)) - (ncol(X) + 1) + 1
  
  tmp <- Beta_extract[ , array_id:(array_id + ncol(X))]                      # Extract single species from list
  
  Beta_posterior[ , j, ] <- t(tmp)          # Fill Beta_posterior
    
}

## R

R_extract <- JSDM$chains$sgibbs                 # Extract R posterior from model (wrong shape)

R_extract <- R_extract[-(1:JSDM$modelList$burnin), ] # Remove burn-in

R_extract <- R_extract[seq(1, nrow(R_extract),       # Implement thinning
                           JSDM$modelList$thin), ]

R_posterior <- array(NA,                        # Create empty array of required shape
                     dim = c(ncol(y),
                             ncol(y),
                             nrow(R_extract)))

for(s in seq_len(n_samples)){                   # Fill correct shape with posterior values
  
  tmp <- matrix(0,                              # 0 not NA so we can use addition later
                nrow = ncol(y),
                ncol = ncol(y))
  
  tmp[lower.tri(tmp,                            # Fill upper.tri of matrix
                diag = TRUE)] <- R_extract[s, ] #  gjam includes diagonal
  
  tmp_diag <- diag(tmp)                         # Save diag because it gets
                                                #  doubled in next step
  tmp <- tmp + t(tmp)                           # Fill lower.tri of matrix
  
  diag(tmp) <- tmp_diag                         # Replace diagonal
  
  tmp <- cov2cor(tmp)                           # Covariance to correlation
  
  R_posterior[ , , s] <- tmp                    # Save each full correlation matrix to an 
                                                #  array slice
}

### Save Posteriors ----

filename <- sprintf("outputs/posteriors/%s_beta_%s_fold_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(Beta_posterior,
        filename)

filename <- sprintf("outputs/posteriors/%s_R_%s_fold_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(R_posterior,
        filename)