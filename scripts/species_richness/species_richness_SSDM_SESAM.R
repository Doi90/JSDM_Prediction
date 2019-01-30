######################################################
######################################################
###                                                ###
###         SSDM SPECIES RICHNESS WORKFLOW         ###
###                                                ###
###    This script calculates the species richness ###
### estimates for the SSDMs                        ###
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

##################################
### Calculate Species Richness ###
##################################

## SSDM - probabilities

### Load Data

filename <- sprintf("outputs/predictions/%s_prob_%s_fold%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

SSDM_pred_prob <- readRDS(filename)

### Calculate Species Richness

SSDM_SR_prob <- species_richness(observed = y_test,
                            predicted = SSDM_pred_prob)

### Save To File

filename <- sprintf("outputs/species_richness/%s_prob_%s_fold%s_SR.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(SSDM_SR_prob,
        filename)

### Memory purge

rm(SSDM_pred_prob,
   SSDM_SR_prob)

## SSDM - binary

### Load Data

filename <- sprintf("outputs/predictions/%s_bin_%s_fold%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

SSDM_pred_bin <- readRDS(filename)

### Calculate Species Richness

SSDM_SR_bin <- species_richness(observed = y_test,
                                 predicted = SSDM_pred_bin)

### Save To File

filename <- sprintf("outputs/species_richness/%s_bin_%s_fold%s_SR.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(SSDM_SR_bin,
        filename)

### Memory purge

rm(SSDM_pred_bin,
   SSDM_SR_bin)

## SESAM

### Load Data

filename <- sprintf("outputs/predictions/SESAM_%s_fold%s.rds",
                    dataset_id,
                    fold_id)

SESAM_pred <- readRDS(filename)

### Calculate Species Richness

SESAM_SR <- species_richness(observed = y_test,
                             predicted = SESAM_pred)

### Save To File

filename <- sprintf("outputs/species_richness/SESAM_%s_fold%s_SR.rds",
                    dataset_id,
                    fold_id)

saveRDS(SESAM_SR,
        filename)

### Memory purge

rm(SESAM_pred,
   SESAM_SR)
