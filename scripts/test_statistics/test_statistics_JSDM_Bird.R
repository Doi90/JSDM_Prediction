######################################################
######################################################
###                                                ###
###         BIRD TEST STATISTICS WORKFLOW          ###
###                                                ###
###    This script calculates the various test     ###
### statistics used in this analysis for the JSDMs ###
### This script is only for the bird runs as it    ###
### splits the jobs up into 10 sample runs         ###
###                                                ###
######################################################
######################################################

#################
### Load Data ###
#################

## Observed / Testing Data

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_test.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y_test <- eval(parse(text = command))                   # Evaluate command to read in data

y_test <- y_test[ , -1]                                 # Remove rownames

#################################
### Calculate Test Statistics ###
#################################

## Marginal ----

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_marginal_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

marg_pred <- readRDS(filename)

### Calculate Test Statistics

marg_ts <- test_statistic(observed = y_test,
                          predictions = marg_pred)

### Save To File

filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_marginal_ts_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(marg_ts,
        filename)

### Memory purge

rm(marg_pred,
   marg_ts)

## Conditional - leave one out ----

### Load data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOO_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

cond_LOO_pred <- readRDS(filename)

### Calculate Test Statistics

cond_LOO_ts <- test_statistic(observed = y_test,
                              predictions = cond_LOO_pred)

### Save To File

filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_condLOO_ts_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(cond_LOO_ts,
        filename)

### Memory purge

rm(cond_LOO_pred,
   cond_LOO_ts)

## Conditional - leave one in ----

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOI_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

cond_LOI_pred <- readRDS(filename)

### Calculate Test Statistics

# Loop over array 4th dimension, calculate test statistic separately
# Save as list, turn back into 4D array

ts_list <- list()

for(i in seq_len(dim(cond_LOI_pred)[4])){
  
  ts_list[[i]] <- test_statistic(observed = y_test,
                                 predictions = cond_LOO_pred[ , , , i])
  
}

cond_LOI_ts <- abind(ts_list,
                     along = 4)

### Save To File

filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_condLOI_ts_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(cond_LOI_ts,
        filename)

### Memory purge

rm(cond_LOI_pred,
   cond_LOI_ts)

## Joint ----

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_joint_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

joint_pred <- readRDS(filename)

### Calculate Test Statistics

joint_ts <- test_statistic(observed = y_test,
                           predictions = joint_pred)

### Save To File

filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_joint_ts_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(joint_ts,
        filename)

### Memory purge

rm(joint_pred,
   joint_ts)

#----