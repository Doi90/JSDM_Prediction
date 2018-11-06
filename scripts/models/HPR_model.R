#########################################################
#########################################################    
###                                                   ###
###                 RUN Pollock JSDM                  ###    
###                                                   ###
### This script runs a Pollock JSDM set up for        ###
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

### Run Pollock JSDM ----

Occur <- as.matrix(y)

X <- as.matrix(X)

## Create regression coefficient matrix

n_species <- ncol(y)

n_env_vars <- ncol(X)

coefs <- matrix(runif(n_species * (n_env_vars + 1), -1, 1),
                ncol=n_env_vars + 1)

## Set up MCMC

n.chains <- 3

n.iter <- 1000000

n.burn <- 15000

n.thin <- 1000

## Set df for inverse-Wishart prior for correlation matrix

df <- 1

## Set name for R object for raw model output

model_name <- "JSDM"

## Run the model

set.seed(28041948)          # Creator's birthday

source('scripts/models/fit_JSDM.R')

### Extract posteriors ----

# Pollock model runs with multiple chains that need to be combined.
# Each subsequent chain is rbind()'d after predecessor.

n_samples <- JSDM[[1]]$BUGSoutput$n.keep        # number of samples per chain

## Beta

Beta_posterior <- array(NA,                     # Create empty array of required shape
                        dim = c(ncol(X) + 1,    # Include intercept
                                ncol(y),
                                n_samples * n.chains))

for(c in seq_len(n.chains)){                               # For each chain
  
  Beta_extract <- JSDM[[c]]$BUGSoutput$sims.list$Beta.raw  # Extract Beta posterior in wrong shape
  
  for(k in seq_len(n_env_vars + 1)){                       # For each covariate including intercept
    
    for(j in seq_len(n_species)){                          # For each species
      
      for(s in seq_len(n_samples)){                        # For each sample in chain
        
        array_dim3_id <- ((c * n_samples) - n_samples) + s  # Need to define a new array
                                                            # Index so successive chains don't overwrite each other
        Beta_posterior[k, j, array_dim3_id] <- Beta_extract[s, j, k]  # Fill empty array with correct values
        
      }
    }
  }
}

## R

R_posterior <- array(NA,                        # Create empty array of required shape
                     dim = c(ncol(y),
                             ncol(y),
                             n_samples * n.chains))

for(c in seq_len(n.chains)){                            # For each chain
  
  R_extract <- JSDM[[c]]$BUGSoutput$sims.list$Tau       # Extract posterior in wrong shape
  
  for(s in seq_len(n_samples)){                         # For each sample
    
    array_dim3_id <- ((c * n_samples) - n_samples) + s  # Need to define a new array
                                                        # Index so successive chains don't overwrite each other
    R_posterior[ , , array_dim3_id] <- cov2cor(R_extract[s, , ])
    
  }
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