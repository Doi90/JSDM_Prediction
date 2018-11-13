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

###########################
### Create Random Folds ###
###########################

# Birds ----

# Bird_y <- Bird[, 1:370]               # Bird PA data
# 
# Bird_X <- Bird[, 371:378]             # Bird covariate data 
# 
# Bird_X <- scale(Bird_X)               # Standardise covariates
# 
# Bird_LL <- Bird[, 379:380]            # Bird coordinate data
# 
# set.seed(28041948)                    # Creator's Birthday
# 
# folds <- createFolds(1:nrow(Bird),                # Number of data points
#                      k=5,                         # 5-fold
#                      list = TRUE,                 # return folds as list
#                      returnTrain = FALSE)         # return all folds, not just train
# 
# for(i in 1:5){
#   
#   fold_id <- paste0("fold",i, collapse = "")                    # set fold id
#   
#   ## Create test/train split for fold i
#   
#   command_y_test <- sprintf('Bird_y[folds$%s,]', fold_id)       # command to call test fold
#   command_y_train <- sprintf('Bird_y[-folds$%s,]', fold_id)     # command to call train folds
#   
#   command_X_test <- sprintf('Bird_X[folds$%s,]', fold_id)       # command to call test fold
#   command_X_train <- sprintf('Bird_X[-folds$%s,]', fold_id)     # command to call train folds
#   
#   command_LL_test <- sprintf('Bird_LL[folds$%s,]', fold_id)     # command to call test fold
#   command_LL_train <- sprintf('Bird_LL[-folds$%s,]', fold_id)   # command to call train folds
#   
#   tmp_y_test <- eval(parse(text = command_y_test))              # evaluate test fold command
#   tmp_y_train <- eval(parse(text = command_y_train))            # evaluate train fold command
#   
#   tmp_X_test <- eval(parse(text = command_X_test))              # evaluate test fold command
#   tmp_X_train <- eval(parse(text = command_X_train))            # evaluate train fold command
#   
#   tmp_LL_test <- eval(parse(text = command_LL_test))            # evaluate test fold command
#   tmp_LL_train <- eval(parse(text = command_LL_train))          # evaluate train fold command
#   
#   ## Save splits to file
#   
#   y_test_filename <- sprintf('data/bird/y_bird_%s_test_random.csv', fold_id)     # test filename
#   y_train_filename <- sprintf('data/bird/y_bird_%s_train_random.csv', fold_id)   # train filename
#   
#   X_test_filename <- sprintf('data/bird/X_bird_%s_test_random.csv', fold_id)     # test filename
#   X_train_filename <- sprintf('data/bird/X_bird_%s_train_random.csv', fold_id)   # train filename
#   
#   LL_test_filename <- sprintf('data/bird/LL_bird_%s_test_random.csv', fold_id)   # test filename
#   LL_train_filename <- sprintf('data/bird/LL_bird_%s_train_random.csv', fold_id) # train filename
#   
#   write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
#   write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
#   
#   write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
#   write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
#   
#   write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
#   write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
#   
# }

# Eucalypts ----

# Euc_y <- Eucalypt[, 1:12]           # Eucalypt PA data
# 
# Euc_X <- Eucalypt[, 13:19]          # Eucalypt covariate data
# 
# Euc_X <- scale(Euc_X)               # Standardise covariates
# 
# Euc_LL <- Eucalypt[, 20:21]         # Eucalypt coordinate data
# 
# set.seed(28041948)                  # Creator's Birthday
# 
# folds <- createFolds(1:nrow(Eucalypt),            # Number of data points
#                      k=5,                         # 5-fold
#                      list = TRUE,                 # return folds as list
#                      returnTrain = FALSE)         # return all folds, not just train
# 
# for(i in 1:5){
#   
#   fold_id <- paste0("fold",i, collapse = "")                    # set fold id
#   
#   command_y_test <- sprintf('Euc_y[folds$%s,]', fold_id)       # command to call test fold
#   command_y_train <- sprintf('Euc_y[-folds$%s,]', fold_id)     # command to call train folds
#   
#   command_X_test <- sprintf('Euc_X[folds$%s,]', fold_id)       # command to call test fold
#   command_X_train <- sprintf('Euc_X[-folds$%s,]', fold_id)     # command to call train folds
#   
#   command_LL_test <- sprintf('Euc_LL[folds$%s,]', fold_id)     # command to call test fold
#   command_LL_train <- sprintf('Euc_LL[-folds$%s,]', fold_id)   # command to call train folds
#   
#   tmp_y_test <- eval(parse(text = command_y_test))              # evaluate test fold command
#   tmp_y_train <- eval(parse(text = command_y_train))            # evaluate train fold command
#   
#   tmp_X_test <- eval(parse(text = command_X_test))              # evaluate test fold command
#   tmp_X_train <- eval(parse(text = command_X_train))            # evaluate train fold command
#   
#   tmp_LL_test <- eval(parse(text = command_LL_test))            # evaluate test fold command
#   tmp_LL_train <- eval(parse(text = command_LL_train))          # evaluate train fold command
#   
#   y_test_filename <- sprintf('data/eucalypt/y_eucalypt_%s_test_random.csv', fold_id)     # test filename
#   y_train_filename <- sprintf('data/eucalypt/y_eucalypt_%s_train_random.csv', fold_id)   # train filename
#   
#   X_test_filename <- sprintf('data/eucalypt/X_eucalypt_%s_test_random.csv', fold_id)     # test filename
#   X_train_filename <- sprintf('data/eucalypt/X_eucalypt_%s_train_random.csv', fold_id)   # train filename
#   
#   LL_test_filename <- sprintf('data/eucalypt/LL_eucalypt_%s_test_random.csv', fold_id)   # test filename
#   LL_train_filename <- sprintf('data/eucalypt/LL_eucalypt_%s_train_random.csv', fold_id) # train filename
#   
#   write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
#   write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
#   
#   write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
#   write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
#   
#   write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
#   write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
#   
# }

# Frogs ----

Frog_y <- Frog[, 1:9]                 # Frog PA data

Frog_X <- Frog[, 10:12]               # Frog covariate data

Frog_X[ , c(1,3)] <- scale(Frog_X[ , c(1,3)]) # Standardise covariates

Frog_LL <- Frog[, 13:14]              # Frog coordinate data

set.seed(28041948)                    # Creator's Birthday

folds <- createFolds(1:nrow(Frog),            # Number of data points
                     k=5,                     # 5-fold
                     list = TRUE,             # return folds as list
                     returnTrain = FALSE)     # return all folds, not just train

for(i in 1:5){
  
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
  
  y_test_filename <- sprintf('data/frog/y_frog_%s_test_spatial.csv', tolower(fold_id))     # test filename
  y_train_filename <- sprintf('data/frog/y_frog_%s_train_spatial.csv', tolower(fold_id))   # train filename
  
  X_test_filename <- sprintf('data/frog/X_frog_%s_test_spatial.csv', tolower(fold_id))     # test filename
  X_train_filename <- sprintf('data/frog/X_frog_%s_train_spatial.csv', tolower(fold_id))   # train filename
  
  LL_test_filename <- sprintf('data/frog/LL_frog_%s_test_spatial.csv', tolower(fold_id))   # test filename
  LL_train_filename <- sprintf('data/frog/LL_frog_%s_train_spatial.csv', tolower(fold_id)) # train filename
  
  write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
  write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
  
  write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
  write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
  
  write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
  write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
  
}

#----

############################
### Create Spatial Folds ###
############################

# Birds ----

Bird_y <- Bird[, 1:370]               # Bird PA data

Bird_X <- Bird[, 371:378]             # Bird covariate data 

Bird_X <- scale(Bird_X)               # Standardise covariates

Bird_LL <- Bird[, 379:380]            # Bird coordinate data

set.seed(28041948)                    # Creator's Birthday

Bird_SP <- SpatialPoints(coords = Bird_LL,
                         proj4string = CRS("+proj=longlat +datum=WGS84"))

Bird_blockCV <- spatialBlock(speciesData = Bird_SP,
                             rows = 5,
                             cols = 10,
                             k = 5)

for(i in 1:5){
  
  fold_id <- paste0("fold",i, collapse = "")          # set fold id
  
  train_id <- Bird_blockCV$folds[[i]][[1]]            # Training site ids
  test_id <- Bird_blockCV$folds[[i]][[2]]             # Testing site ids
  
  tmp_y_test <- Bird_y[test_id, ]                     # y test data
  tmp_y_train <- Bird_y[train_id, ]                   # y train data
  
  tmp_X_test <- Bird_X[test_id, ]                     # X test data
  tmp_X_train <- Bird_X[train_id, ]                   # X train data
  
  tmp_LL_test <- Bird_LL[test_id, ]                   # Lat/long test data
  tmp_LL_train <- Bird_LL[train_id, ]                 # Lat/long train data
  
  y_test_filename <- sprintf('data/bird/y_bird_%s_test_spatial.csv', fold_id)     # test filename
  y_train_filename <- sprintf('data/bird/y_bird_%s_train_spatial.csv', fold_id)   # train filename
  
  X_test_filename <- sprintf('data/bird/X_bird_%s_test_spatial.csv', fold_id)     # test filename
  X_train_filename <- sprintf('data/bird/X_bird_%s_train_spatial.csv', fold_id)   # train filename
  
  LL_test_filename <- sprintf('data/bird/LL_bird_%s_test_spatial.csv', fold_id)   # test filename
  LL_train_filename <- sprintf('data/bird/LL_bird_%s_train_spatial.csv', fold_id) # train filename
  
  write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
  write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
  
  write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
  write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
  
  write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
  write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
  
}

# Eucalypts ----

Euc_y <- Eucalypt[, 1:12]           # Eucalypt PA data

Euc_X <- Eucalypt[, 13:19]          # Eucalypt covariate data

Euc_X <- scale(Euc_X)               # Standardise covariates

Euc_LL <- Eucalypt[, 20:21]         # Eucalypt coordinate data

set.seed(28041948)                  # Creator's Birthday

Euc_SP <- SpatialPoints(coords = Euc_LL,
                        proj4string = CRS("+proj=longlat +datum=WGS84"))

Euc_blockCV <- spatialBlock(speciesData = Euc_SP,
                            rows = 10,
                            cols = 5,
                            k = 5)

for(i in 1:5){
  
  fold_id <- paste0("fold",i, collapse = "")          # set fold id
  
  train_id <- Euc_blockCV$folds[[i]][[1]]             # Training site ids
  test_id <- Euc_blockCV$folds[[i]][[2]]              # Testing site ids
  
  tmp_y_test <- Euc_y[test_id, ]                      # y test data
  tmp_y_train <- Euc_y[train_id, ]                    # y train data
  
  tmp_X_test <- Euc_X[test_id, ]                      # X test data
  tmp_X_train <- Euc_X[train_id, ]                    # X train data
  
  tmp_LL_test <- Euc_LL[test_id, ]                    # Lat/long test data
  tmp_LL_train <- Euc_LL[train_id, ]                  # Lat/long train data
  
  y_test_filename <- sprintf('data/eucalypt/y_eucalypt_%s_test_spatial.csv', fold_id)     # test filename
  y_train_filename <- sprintf('data/eucalypt/y_eucalypt_%s_train_spatial.csv', fold_id)   # train filename
  
  X_test_filename <- sprintf('data/eucalypt/X_eucalypt_%s_test_spatial.csv', fold_id)     # test filename
  X_train_filename <- sprintf('data/eucalypt/X_eucalypt_%s_train_spatial.csv', fold_id)   # train filename
  
  LL_test_filename <- sprintf('data/eucalypt/LL_eucalypt_%s_test_spatial.csv', fold_id)   # test filename
  LL_train_filename <- sprintf('data/eucalypt/LL_eucalypt_%s_train_spatial.csv', fold_id) # train filename
  
  write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
  write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
  
  write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
  write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
  
  write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
  write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
  
}

# Frogs ----

# Frog_y <- Frog[, 1:9]                 # Frog PA data
# 
# Frog_X <- Frog[, 10:12]               # Frog covariate data
# 
# Frog_X <- scale(Frog_X)               # Standardise covariates
# 
# Frog_LL <- Frog[, 13:14]              # Frog coordinate data
# 
# set.seed(28041948)                    # Creator's Birthday
# 
# Frog_SP <- SpatialPoints(coords = Frog_LL,
#                          proj4string = CRS("+proj=longlat +datum=WGS84"))
# 
# Frog_blockCV <- spatialBlock(speciesData = Frog_SP,
#                              rows = 1,
#                              cols = 2,
#                              k = 2)
# for(i in 1:2){
#   
#   fold_id <- paste0("fold",i, collapse = "")          # set fold id
#   
#   train_id <- Frog_blockCV$folds[[i]][[1]]            # Training site ids
#   test_id <- Frog_blockCV$folds[[i]][[2]]             # Testing site ids
#   
#   tmp_y_test <- Frog_y[test_id, ]                     # y test data
#   tmp_y_train <- Frog_y[train_id, ]                   # y train data
#   
#   tmp_X_test <- Frog_X[test_id, ]                     # X test data
#   tmp_X_train <- Frog_X[train_id, ]                   # X train data
#   
#   tmp_LL_test <- Frog_LL[test_id, ]                   # Lat/long test data
#   tmp_LL_train <- Frog_LL[train_id, ]                 # Lat/long train data
#   
#   y_test_filename <- sprintf('data/frog/y_frog_%s_test_spatial.csv', fold_id)     # test filename
#   y_train_filename <- sprintf('data/frog/y_frog_%s_train_spatial.csv', fold_id)   # train filename
#   
#   X_test_filename <- sprintf('data/frog/X_frog_%s_test_spatial.csv', fold_id)     # test filename
#   X_train_filename <- sprintf('data/frog/X_frog_%s_train_spatial.csv', fold_id)   # train filename
#   
#   LL_test_filename <- sprintf('data/frog/LL_frog_%s_test_spatial.csv', fold_id)   # test filename
#   LL_train_filename <- sprintf('data/frog/LL_frog_%s_train_spatial.csv', fold_id) # train filename
#   
#   write.csv(tmp_y_test, y_test_filename)              # write y fold i test csv
#   write.csv(tmp_y_train, y_train_filename)            # write y fold i train csv
#   
#   write.csv(tmp_X_test, X_test_filename)              # write X fold i test csv
#   write.csv(tmp_X_train, X_train_filename)            # write X fold i train csv
#   
#   write.csv(tmp_LL_test, LL_test_filename)            # write LL fold i test csv
#   write.csv(tmp_LL_train, LL_train_filename)          # write LL fold i train csv
#   
# }
