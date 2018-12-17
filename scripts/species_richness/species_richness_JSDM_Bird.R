######################################################
######################################################
###                                                ###
###         BIRD SPECIES RICHNESS WORKFLOW         ###
###                                                ###
###    This script calculates the species richness ###
### estimates for the JSDMs. This script is only   ###
### for the bird runs as it runs on only a subset  ###
### of posterior samples (10).                     ###
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

## Marginal 

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_marginal_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

marg_pred <- readRDS(filename)

### Calculate species richness

marg_SR <- species_richness(observed = y_test,
                            predicted = marg_pred)

### Save To File

filename <- sprintf("outputs/species_richness/%s_%s_fold%s_marginal_SR_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(marg_SR,
        filename)

### Memory purge

rm(marg_pred,
   marg_SR)

## Conditional - leave one out

### Load data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOO_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

cond_LOO_pred <- readRDS(filename)

### Calculate species richness

cond_LOO_SR <- species_richness(observed = y_test,
                                predicted = cond_LOO_pred)

### Save To File

filename <- sprintf("outputs/species_richness/%s_%s_fold%s_condLOO_SR_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(cond_LOO_SR,
        filename)

### Memory purge

rm(cond_LOO_pred,
   cond_LOO_SR)

## Conditional - leave one in

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOI_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

cond_LOI_pred <- readRDS(filename)

### Calculate species richness

SR_list <- list()

for(i in seq_len(dim(cond_LOI_pred)[4])){
  
  SR_list[[i]] <- species_richness(observed = y_test,
                                   predicted = cond_LOI_pred[ , , , i])
  
}

cond_LOI_SR <- abind(SR_list,
                     along = 4)

### Save To File

filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_condLOI_SR_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(cond_LOI_SR,
        filename)

### Memory purge

rm(cond_LOI_pred,
   cond_LOI_tSR)

## Joint

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_joint_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

joint_pred <- readRDS(filename)

### Calculate species richness

joint_SR <- species_richness(observed = y_test,
                             predicted = joint_pred)

### Save To File

filename <- sprintf("outputs/species_richness/%s_%s_fold%s_joint_SR_%s_%s.rds",
                    model_id,
                    dataset_id,
                    fold_id,
                    start_sample,
                    end_sample)

saveRDS(joint_SR,
        filename)

### Memory purge

rm(joint_pred,
   joint_SR)

