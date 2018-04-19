#########################################################
#########################################################    
###                                                   ###
###               Stacked probit GLM                  ###    
###                                                   ###
### This script runs a Stacked probit GLM set up for  ###
### cross-validation on Spartan. Also defines outputs ###
### in correct format for prediction code.            ###  
###                                                   ###
#########################################################
#########################################################

### Load Packages ----

library(coda)

### Load Data ----

## CV fold id

fold_id <- read.table("fold.txt")[1,1]   # CV fold id needs to be read in from text file.
                                         # Pass in Spartan index.
## Dataset id

dataset_id <- read.table("dataset.txt",                 # Dataset id needs to be read in from
                         stringsAsFactors = FALSE)[1,1] # text file

## Pres/Abs data

command <- sprintf("read.csv('y_%s_Fold%s_train.csv')", # Need to build command to read in
                   dataset_id,                          # specific files for this CV fold
                   fold_id)

y <- eval(parse(text = command))                        # Evaluate command to read in data

y <- y[ , -1]                                           # Remove rownames

## Site data

command <- sprintf("read.csv('x_%s_Fold%s_train.csv')", # Need to build command to read in
                   dataset_id,                          # specific files for this CV fold
                   fold_id)

X <- eval(parse(text = command))                        # Evaluate command to read in data

X <- X[ , -1]                                           # Remove rownames

### Run Stacked probit GLM ----

set.seed(28041948)                                      # Creator's Birthday

model_list <- list()                                    # Empty list to store models

for(j in seq_len(ncol(y))){                             # For each species
  
  data <- cbind(y[ , j],                                # Create single species dataframe
                X)
  
  colnames(data)[1] <- "Occurrence"                     # Standardise P/A column name
  
  form <- as.formula(paste("Occurrence ~ ",             # Build formula from X column names
                           paste(colnames(X)[-1],
                                 collapse = "+")))
  
  ## Build command to run glm and save to list
  
  command <- sprintf("model_list[[j]] <-                 
                         glm_%s <- glm(formula = form,
                                       data = data,
                                       family = binomial(link = 'probit')",
                     colnames(y)[j])
  
  eval(parse(text = command))                           # Run command
  
}

## Save Model Outputs ----

filename <- sprintf("SSDM_%s_fold_%s.rds",
                    dataset_id,
                    fold_id)
saveRDS(model_list,
        filename)