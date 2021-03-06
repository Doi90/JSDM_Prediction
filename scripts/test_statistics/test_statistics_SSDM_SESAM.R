######################################################
######################################################
###                                                ###
###         SSDM TEST STATISTICS WORKFLOW          ###
###                                                ###
###    This script calculates the various test     ###
### statistics used in this analysis for the SSDMs ###
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

## SSDM - probabilities

### Load Data

filename <- sprintf("outputs/predictions/%s_prob_%s_fold%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

SSDM_pred_prob <- readRDS(filename)

### Calculate Test Statistics

SSDM_ts_prob <- test_statistic(observed = y_test,
                          predictions = SSDM_pred_prob)

### Save To File

filename <- sprintf("outputs/test_statistics/%s_prob_%s_fold%s_ts.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(SSDM_ts_prob,
        filename)

### Memory purge

rm(SSDM_pred_prob,
   SSDM_ts_prob)

## SSDM - binary

### Load Data

filename <- sprintf("outputs/predictions/%s_bin_%s_fold%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

SSDM_pred_bin <- readRDS(filename)

### Calculate Test Statistics

SSDM_ts_bin <- test_statistic(observed = y_test,
                               predictions = SSDM_pred_bin)

### Save To File

filename <- sprintf("outputs/test_statistics/%s_bin_%s_fold%s_ts.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(SSDM_ts_bin,
        filename)

### Memory purge

rm(SSDM_pred_bin,
   SSDM_ts_bin)

## SESAM

### Load Data

filename <- sprintf("outputs/predictions/SESAM_%s_fold%s.rds",
                    dataset_id,
                    fold_id)

SESAM_pred <- readRDS(filename)

### Calculate Test Statistics

SESAM_ts <- test_statistic(observed = y_test,
                           predictions = SESAM_pred)

### Save To File

filename <- sprintf("outputs/test_statistics/SESAM_%s_fold%s_ts.rds",
                    dataset_id,
                    fold_id)

saveRDS(SESAM_ts,
        filename)

### Memory purge

rm(SESAM_pred,
   SESAM_ts)
