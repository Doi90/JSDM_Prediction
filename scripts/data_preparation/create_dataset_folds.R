######################################################
######################################################
###                                                ###
###              Create Dataset Folds              ###
###                                                ###
###  This script generates the test/train splits   ###
### for random and spatial cross-validation folds. ###
### Needs to be done manually before model fitting ###
### to ensure they are identical between model     ###
### runs on Spartan.                               ###
###                                                ###
######################################################
######################################################

#####################
### Load Packages ###
#####################

library(caret)
library(blockCV)
library(sp)

#################
### Load Data ###
#################

Bird <- read.csv("data/bird/birds_compiled.csv")
Eucalypt <- read.csv("data/eucalypt/eucalypt_compiled.csv")
Frog <- read.csv("data/frog/frog_compiled.csv")
Butterfly <- read.csv("data/butterfly/Butterfly_Compiled.csv")

###############################
### Create Custom Functions ###
###############################

## colMax() gets the maximum value in all columns

colMax <- function(df){
  
  max_values <- vector(length = ncol(df))
  
  for(i in seq_len(ncol(df))){
    
    max_values[i] <- max(df[ , i])
    
  }
  
  return(max_values)
  
}

## blockCV_multi_spp_test() checks the generated folds
## meet our requirements. Training data in each CV split
## must have at least one presence for all species.

blockCV_multi_spp_test <- function(dataset,
                                   blockCV){
  
  train_tests <- vector(length = 5)
  
  for(i in seq_len(5)){
    
    train_id <- blockCV$folds[[i]][[1]]
    
    y_train <- dataset[train_id, ]
    
    train_tests[i] <- all(colMax(y_train) > 0)
    
  }
  
  if(all(train_tests == TRUE)){
    
    return(TRUE)
    
  } else {
    
    return(FALSE)
  }
  
}

## blockCV_multi_spp_test() checks the generated folds
## meet our requirements. Training data in each CV split
## must have at least one presence for all species.

caret_multi_spp_test <- function(dataset,
                                 folds){
  
  train_tests <- vector(length = 5)
  
  for(i in seq_len(5)){
    
    fold_id <- paste0("Fold", i, collapse = "")
    
    command <- sprintf("folds$%s", fold_id)
    
    train_id <- eval(parse(text = command)) # rows in TEST DATA
    
    y_train <- dataset[-train_id, ]
    
    train_tests[i] <- all(colMax(y_train) > 0)
    
  }
  
  if(all(train_tests == TRUE)){
    
    return(TRUE)
    
  } else {
    
    return(FALSE)
  }
  
}

######################
### Generate Folds ###
######################

# Birds ----

## Data manipulation

Bird_y <- Bird[, 1:370]                           # Bird PA data

Bird_y <- Bird_y[, -which(colSums(Bird_y) < 50)]  # Only keep species with 50+ presences

Bird_X <- Bird[, 371:378]                         # Bird covariate data 

Bird_X <- scale(Bird_X)                           # Standardise covariates

Bird_LL <- Bird[, 379:380]                        # Bird coordinate data

## Generate folds

### Define control parameters

dataset_ok <- FALSE

iteration <- 0

### Generate folds until meet our test condition

while(dataset_ok == FALSE){
  
  ## Increase counter
  
  iteration <- iteration + 1
  
  ## Prevent infinite loop
  
  if(iteration >= 1000){
    
    break()
    
  }
  
  ## Generate folds
  
  set.seed(28041948 + iteration)                    # Creator's Birthday
  
  Bird_SP <- SpatialPoints(coords = Bird_LL,
                           proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  Bird_blockCV <- spatialBlock(speciesData = Bird_SP,
                               rows = 20,
                               cols = 40,
                               k = 5,
                               iteration = 100)
  
  dataset_ok <- blockCV_multi_spp_test(dataset = Bird_y,
                                       blockCV = Bird_blockCV)
  
}

## Write data to file

bird_site_ids <- list()

for(i in seq_len(5)){
  
  fold_id <- paste0("fold",i, collapse = "")          # set fold id
  
  train_id <- Bird_blockCV$folds[[i]][[1]]            # Training site ids
  test_id <- Bird_blockCV$folds[[i]][[2]]             # Testing site ids
  
  tmp_y_test <- Bird_y[test_id, ]                     # y test data
  tmp_y_train <- Bird_y[train_id, ]                   # y train data
  
  tmp_X_test <- Bird_X[test_id, ]                     # X test data
  tmp_X_train <- Bird_X[train_id, ]                   # X train data
  
  tmp_LL_test <- Bird_LL[test_id, ]                   # Lat/long test data
  tmp_LL_train <- Bird_LL[train_id, ]                 # Lat/long train data
  
  y_test_filename <- sprintf('data/bird/y_bird_%s_test.csv', fold_id)     # test filename
  y_train_filename <- sprintf('data/bird/y_bird_%s_train.csv', fold_id)   # train filename
  
  X_test_filename <- sprintf('data/bird/X_bird_%s_test.csv', fold_id)     # test filename
  X_train_filename <- sprintf('data/bird/X_bird_%s_train.csv', fold_id)   # train filename
  
  LL_test_filename <- sprintf('data/bird/LL_bird_%s_test.csv', fold_id)   # test filename
  LL_train_filename <- sprintf('data/bird/LL_bird_%s_train.csv', fold_id) # train filename
  
  write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
  write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
  
  write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
  write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
  
  write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
  write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
  
  bird_site_ids[[i]] <- Bird_blockCV$folds[[i]][[2]]
  
}

saveRDS(bird_site_ids,
        "data/bird/site_ids.rds")

# Eucalypts ----

## Data manipulation

Euc_y <- Eucalypt[, 1:12]           # Eucalypt PA data

Euc_X <- Eucalypt[, 13:19]          # Eucalypt covariate data

Euc_X[ , -c(2, 5)] <- scale(Euc_X[ , -c(2, 5)])  # Standardise covariates

Euc_LL <- Eucalypt[, 20:21]         # Eucalypt coordinate data

## Generate folds

### Define control parameters

dataset_ok <- FALSE

iteration <- 0

### Generate folds until meet our test condition

while(dataset_ok == FALSE){
  
  ## Increase counter
  
  iteration <- iteration + 1
  
  ## Prevent infinite loop
  
  if(iteration >= 1000){
    
    break()
    
  }
  
  ## Generate folds
  
  set.seed(28041948 + iteration)                    # Creator's Birthday
  
  Euc_SP <- SpatialPoints(coords = Euc_LL,
                           proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  Euc_blockCV <- spatialBlock(speciesData = Euc_SP,
                               rows = 10,
                               cols = 5,
                               k = 5,
                               iteration = 100)
  
  dataset_ok <- blockCV_multi_spp_test(dataset = Euc_y,
                                       blockCV = Euc_blockCV)
  
}

## Write data to file

euc_site_ids <- list()

for(i in seq_len(5)){
  
  fold_id <- paste0("fold",i, collapse = "")          # set fold id
  
  train_id <- Euc_blockCV$folds[[i]][[1]]             # Training site ids
  test_id <- Euc_blockCV$folds[[i]][[2]]              # Testing site ids
  
  tmp_y_test <- Euc_y[test_id, ]                      # y test data
  tmp_y_train <- Euc_y[train_id, ]                    # y train data
  
  tmp_X_test <- Euc_X[test_id, ]                      # X test data
  tmp_X_train <- Euc_X[train_id, ]                    # X train data
  
  tmp_LL_test <- Euc_LL[test_id, ]                    # Lat/long test data
  tmp_LL_train <- Euc_LL[train_id, ]                  # Lat/long train data
  
  y_test_filename <- sprintf('data/eucalypt/y_eucalypt_%s_test.csv', fold_id)     # test filename
  y_train_filename <- sprintf('data/eucalypt/y_eucalypt_%s_train.csv', fold_id)   # train filename
  
  X_test_filename <- sprintf('data/eucalypt/X_eucalypt_%s_test.csv', fold_id)     # test filename
  X_train_filename <- sprintf('data/eucalypt/X_eucalypt_%s_train.csv', fold_id)   # train filename
  
  LL_test_filename <- sprintf('data/eucalypt/LL_eucalypt_%s_test.csv', fold_id)   # test filename
  LL_train_filename <- sprintf('data/eucalypt/LL_eucalypt_%s_train.csv', fold_id) # train filename
  
  write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
  write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
  
  write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
  write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
  
  write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
  write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
  
  euc_site_ids[[i]] <- Euc_blockCV$folds[[i]][[2]]
  
}

saveRDS("data/eucalypt/site_ids.rds")

# Frogs ----

## Data manipulation

Frog_y <- Frog[, 1:9]                 # Frog PA data

Frog_X <- Frog[, 10:12]               # Frog covariate data

Frog_X[ , c(1,3)] <- scale(Frog_X[ , c(1,3)]) # Standardise covariates

Frog_LL <- Frog[, 13:14]              # Frog coordinate data

## Generate folds

### Define control parameters

dataset_ok <- FALSE

iteration <- 0

### Generate folds until meet our test condition

while(dataset_ok == FALSE){
  
  ## Increase counter
  
  iteration <- iteration + 1
  
  ## Prevent infinite loop
  
  if(iteration >= 1000){
    
    break()
    
  }
  
  ## Generate folds
  
  set.seed(28041948 + iteration)                    # Creator's Birthday
  
  folds <- createFolds(1:nrow(Frog),            # Number of data points
                       k=5,                     # 5-fold
                       list = TRUE,             # return folds as list
                       returnTrain = FALSE)     # return all folds, not just train
  
  dataset_ok <- caret_multi_spp_test(dataset = Frog_y,
                                     folds = folds)
  
}

frog_site_ids <- list()

for(i in seq_len(5)){
  
  fold_id <- paste0("Fold",i, collapse = "")                    # set fold id
  
  command_y_test <- sprintf('Frog_y[folds$%s,]', fold_id)       # command to call test fold
  command_y_train <- sprintf('Frog_y[-folds$%s,]', fold_id)     # command to call train folds
  
  command_X_test <- sprintf('Frog_X[folds$%s,]', fold_id)       # command to call test fold
  command_X_train <- sprintf('Frog_X[-folds$%s,]', fold_id)     # command to call train folds
  
  command_LL_test <- sprintf('Frog_LL[folds$%s,]', fold_id)     # command to call test fold
  command_LL_train <- sprintf('Frog_LL[-folds$%s,]', fold_id)   # command to call train folds
  
  tmp_y_test <- eval(parse(text = command_y_test))              # evaluate test fold command
  tmp_y_train <- eval(parse(text = command_y_train))            # evaluate train fold command
  
  tmp_X_test <- eval(parse(text = command_X_test))              # evaluate test fold command
  tmp_X_train <- eval(parse(text = command_X_train))            # evaluate train fold command
  
  tmp_LL_test <- eval(parse(text = command_LL_test))            # evaluate test fold command
  tmp_LL_train <- eval(parse(text = command_LL_train))          # evaluate train fold command
  
  y_test_filename <- sprintf('data/frog/y_frog_%s_test.csv', tolower(fold_id))     # test filename
  y_train_filename <- sprintf('data/frog/y_frog_%s_train.csv', tolower(fold_id))   # train filename
  
  X_test_filename <- sprintf('data/frog/X_frog_%s_test.csv', tolower(fold_id))     # test filename
  X_train_filename <- sprintf('data/frog/X_frog_%s_train.csv', tolower(fold_id))   # train filename
  
  LL_test_filename <- sprintf('data/frog/LL_frog_%s_test.csv', tolower(fold_id))   # test filename
  LL_train_filename <- sprintf('data/frog/LL_frog_%s_train.csv', tolower(fold_id)) # train filename
  
  write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
  write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
  
  write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
  write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
  
  write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
  write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
  
  frog_site_ids[[i]] <- folds[[i]]
  
}

saveRDS(frog_site_ids,
        "data/frog/site_ids.rds")

#----

# Butterflies ----

## Data manipulation

Butterfly_y <- Butterfly[ , 1:55]                 # Butterfly PA data

Butterfly_y <- Butterfly_y[ , -which(colSums(Butterfly_y) < 50)]  # Only keep species with 50+ presences

Butterfly_X <- Butterfly[ , 56:59]                # Butterfly covariate data 

Butterfly_X <- scale(Butterfly_X)                 # Standardise covariates

Butterfly_LL <- Butterfly[, 60:61]                # Butterfly coordinate data

## Generate folds

### Define control parameters

dataset_ok <- FALSE

iteration <- 0

### Generate folds until meet our test condition

while(dataset_ok == FALSE){
  
  ## Increase counter
  
  iteration <- iteration + 1
  
  ## Prevent infinite loop
  
  if(iteration >= 1000){
    
    break()
    
  }
  
  ## Generate folds
  
  set.seed(28041948 + iteration)                    # Creator's Birthday
  
  Butterfly_SP <- SpatialPoints(coords = Butterfly_LL,
                                proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  Butterfly_blockCV <- spatialBlock(speciesData = Butterfly_SP,
                                    rows = 20,
                                    cols = 10,
                                    k = 5,
                                    iteration = 100)
  
  dataset_ok <- blockCV_multi_spp_test(dataset = Butterfly_y,
                                       blockCV = Butterfly_blockCV)
  
}

## Write data to file

butterfly_site_ids <- list()

for(i in seq_len(5)){
  
  fold_id <- paste0("fold",i, collapse = "")          # set fold id
  
  train_id <- Butterfly_blockCV$folds[[i]][[1]]            # Training site ids
  test_id <- Butterfly_blockCV$folds[[i]][[2]]             # Testing site ids
  
  tmp_y_test <- Butterfly_y[test_id, ]                     # y test data
  tmp_y_train <- Butterfly_y[train_id, ]                   # y train data
  
  tmp_X_test <- Butterfly_X[test_id, ]                     # X test data
  tmp_X_train <- Butterfly_X[train_id, ]                   # X train data
  
  tmp_LL_test <- Butterfly_LL[test_id, ]                   # Lat/long test data
  tmp_LL_train <- Butterfly_LL[train_id, ]                 # Lat/long train data
  
  y_test_filename <- sprintf('data/butterfly/y_butterfly_%s_test.csv', fold_id)     # test filename
  y_train_filename <- sprintf('data/butterfly/y_butterfly_%s_train.csv', fold_id)   # train filename
  
  X_test_filename <- sprintf('data/butterfly/X_butterfly_%s_test.csv', fold_id)     # test filename
  X_train_filename <- sprintf('data/butterfly/X_butterfly_%s_train.csv', fold_id)   # train filename
  
  LL_test_filename <- sprintf('data/butterfly/LL_butterfly_%s_test.csv', fold_id)   # test filename
  LL_train_filename <- sprintf('data/butterfly/LL_butterfly_%s_train.csv', fold_id) # train filename
  
  write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
  write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
  
  write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
  write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
  
  write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
  write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
  
  butterfly_site_ids[[i]] <- Butterfly_blockCV$folds[[i]][[2]]
  
}

saveRDS(butterfly_site_ids,
        "data/butterfly/site_ids.rds")
