#########################################################
#########################################################    
###                                                   ###
###               STACKED PROBIT GLM                  ###    
###                                                   ###
### This script runs a Stacked probit GLM set up for  ###
### cross-validation on Spartan. Also defines outputs ###
### in correct format for prediction code.            ###  
###                                                   ###
#########################################################
#########################################################

## Load Data ----

### Pres/Abs data

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_train_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y <- eval(parse(text = command))                        # Evaluate command to read in data

y <- y[ , -1]                                           # Remove rownames

### Site data

command <- sprintf("read.csv('data/%1$s/X_%1$s_fold%2$s_train_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

X <- eval(parse(text = command))                        # Evaluate command to read in data

X <- X[ , -1]                                           # Remove rownames

## Run Stacked probit GLM ----

set.seed(28041948)                                      # Creator's Birthday

model_list <- vector("list", ncol(y))                   # Empty list to store models

for(j in seq_len(ncol(y))){                             # For each species
  
  data <- cbind(y[ , j],                                # Create single species dataframe
                X)
  
  colnames(data)[1] <- "Occurrence"                     # Standardise P/A column name
  
  form <- as.formula(paste("Occurrence ~ ",             # Build formula from X column names
                           paste(colnames(X)[-1],
                                 collapse = "+")))
  
  ## Build command to run glm and save to list
  
  command <- sprintf("model_list[[%s]] <-                 
                         glm_%s <- glm(formula = form,
                                       data = data,
                                       family = binomial(link = 'probit')",
                     j,
                     colnames(y)[j])
  
  eval(parse(text = command))                           # Run command
  
}

## Save Model Outputs ----

filename <- sprintf("posteriors/%s_%s_fold_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(model_list,
        filename)
