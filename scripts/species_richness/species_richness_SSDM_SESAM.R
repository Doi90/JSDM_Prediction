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

command <- sprintf("read.csv('data/%1$s/y_%1$s_fold%2$s_test_spatial.csv')", 
                   dataset_id,                          # Need to build command to read in
                   fold_id)                             # specific files for this CV fold

y_test <- eval(parse(text = command))                   # Evaluate command to read in data

y_test <- y_test[ , -1]                                 # Remove rownames

##################################
### Calculate Species Richness ###
##################################

## SSDM

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s.rds",
                    model_id,
                    dataset_id,
                    fold_id)

SSDM_pred <- readRDS(filename)

### Calculate Species Richness

SSDM_SR <- species_richness(observed = y_test,
                            predictions = SSDM_pred)

### Save To File

filename <- sprintf("outputs/species_richness/%s_%s_fold%s_SR.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(SSDM_SR,
        filename)

### Memory purge

rm(SSDM_pred,
   SSDM_SR)

## SESAM

### Load Data

filename <- sprintf("outputs/predictions/SESAM_%s_fold%s.rds",
                    dataset_id,
                    fold_id)

SESAM_pred <- readRDS(filename)

### Calculate Species Richness

SESAM_SR <- species_richness(observed = y_test,
                             predictions = SESAM_pred)

### Save To File

filename <- sprintf("outputs/species_richness/SESAM_%s_fold%s_SR.rds",
                    dataset_id,
                    fold_id)

saveRDS(SESAM_SR,
        filename)

### Memory purge

rm(SESAM_pred,
   SESAM_SR)
