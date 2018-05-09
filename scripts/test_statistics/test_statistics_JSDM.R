######################################################
######################################################
###                                                ###
###           Test Statistics Workflow             ###
###                                                ###
###    This script calculates the various test     ###
###       statistics used in this analysis.        ###
###                                                ###
######################################################
######################################################

#################
### Load Data ###
#################

## Observed / Testing Data

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_test_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y_test <- eval(parse(text = command))                   # Evaluate command to read in data

y_test <- y_test[ , -1]                                 # Remove rownames

#################################
### Calculate Test Statistics ###
#################################

## Marginal 

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_marginal.rds",
                    model_id,
                    dataset_id,
                    fold_id)

marg_pred <- readRDS(filename)

### Calculate Test Statistics

marg_ts <- test_statistic(observed = y_test,
                          predictions = marg_pred)

### Save To File

filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_marginal_ts.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(marg_ts,
        filename)

### Memory purge

rm(c("marg_pred",
     "marg_ts"))

## Conditional - leave one out

### Load data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOO.rds",
                    model_id,
                    dataset_id,
                    fold_id)

cond_LOO_pred <- readRDS(filename)

### Calculate Test Statistics

cond_LOO_ts <- test_statistic(observed = y_test,
                              predictions = cond_LOO_pred)

### Save To File

filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_cond_LOO_ts.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(cond_LOO_ts,
        filename)

### Memory purge

rm(c("cond_LOO_pred",
     "cond_LOO_ts"))

## Conditional - leave one in

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOI.rds",
                    model_id,
                    dataset_id,
                    fold_id)

cond_LOI_pred <- readRDS(filename)

### Calculate Test Statistics

# Loop over array 4th dimension, fill new empty array, then average over that?

### Save To File

filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_cond_LOI_ts.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(cond_LOI_ts,
        filename)

### Memory purge

rm(c("cond_LOI_pred",
     "cond_LOI_ts"))

## Joint

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_joint.rds",
                    model_id,
                    dataset_id,
                    fold_id)

joint_pred <- readRDS(filename)

### Calculate Test Statistics

joint_ts <- test_statistic(observed = y_test,
                           predictions = joint_pred)

### Save To File

filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_joint_ts.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(joint_ts,
        filename)

### Memory purge

rm(c("joint_pred",
     "joint_ts"))
