######################################################
######################################################
###                                                ###
###         JSDM SPECIES RICHNESS WORKFLOW         ###
###                                                ###
###    This script calculates the species richness ###
### estimates for the JSDMs                        ###
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

filename <- sprintf("outputs/predictions/%s_%s_fold%s_marginal.rds",
                    model_id,
                    dataset_id,
                    fold_id)

marg_pred <- readRDS(filename)

### Calculate species richness

marg_SR <- species_richness(observed = y_test,
                            predicted = marg_pred)

### Save To File

filename <- sprintf("outputs/species_richness/%s_%s_fold%s_marginal_SR.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(marg_SR,
        filename)

### Memory purge

rm(marg_pred,
   marg_SR)

## Conditional - leave one out

### Load data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOO.rds",
                    model_id,
                    dataset_id,
                    fold_id)

cond_LOO_pred <- readRDS(filename)

### Calculate species richness

cond_LOO_SR <- species_richness(observed = y_test,
                                predicted = cond_LOO_pred)

### Save To File

filename <- sprintf("outputs/species_richness/%s_%s_fold%s_condLOO_SR.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(cond_LOO_SR,
        filename)

### Memory purge

rm(cond_LOO_pred,
   cond_LOO_SR)

# ## Conditional - leave one in
# 
# ### Load Data
# 
# filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOI.rds",
#                     model_id,
#                     dataset_id,
#                     fold_id)
# 
# cond_LOI_pred <- readRDS(filename)
# 
# ### Calculate species richness
# 
# # Loop over array 4th dimension, fill new empty array, then average over that?
# 
# ### Save To File
# 
# filename <- sprintf("outputs/test_statistics/%s_%s_fold%s_condLOI_ts.rds",
#                     model_id,
#                     dataset_id,
#                     fold_id)
# 
# saveRDS(cond_LOI_ts,
#         filename)
# 
# ### Memory purge
# 
# rm(cond_LOI_pred,
#    cond_LOI_ts)

## Joint

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_joint.rds",
                    model_id,
                    dataset_id,
                    fold_id)

joint_pred <- readRDS(filename)

### Calculate species richness

joint_SR <- species_richness(observed = y_test,
                             predicted = joint_pred)

### Save To File

filename <- sprintf("outputs/species_richness/%s_%s_fold%s_joint_SR.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(joint_SR,
        filename)

### Memory purge

rm(joint_pred,
   joint_SR)

