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

## Marginal - probabilities

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_marginal_prob.rds",
                    model_id,
                    dataset_id,
                    fold_id)

marg_pred_prob <- readRDS(filename)

### Calculate species richness

marg_SR_prob <- species_richness(observed = y_test,
                                 predicted = marg_pred_prob)

### Save To File

filename <- sprintf("outputs/species_richness/%s_%s_fold%s_marginal_prob_SR.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(marg_SR_prob,
        filename)

### Memory purge

rm(marg_pred_prob,
   marg_SR_prob)

## Marginal - binary

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_marginal_bin.rds",
                    model_id,
                    dataset_id,
                    fold_id)

marg_pred_bin <- readRDS(filename)

### Calculate species richness

marg_SR_bin <- species_richness(observed = y_test,
                                predicted = marg_pred_bin)

### Save To File

filename <- sprintf("outputs/species_richness/%s_%s_fold%s_marginal_bin_SR.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(marg_SR_bin,
        filename)

### Memory purge

rm(marg_pred_bin,
   marg_SR_bin)

# ## Conditional - leave one out
# 
# ### Load data
# 
# filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOO.rds",
#                     model_id,
#                     dataset_id,
#                     fold_id)
# 
# cond_LOO_pred <- readRDS(filename)
# 
# ### Calculate species richness
# 
# cond_LOO_SR <- species_richness(observed = y_test,
#                                 predicted = cond_LOO_pred)
# 
# ### Save To File
# 
# filename <- sprintf("outputs/species_richness/%s_%s_fold%s_condLOO_SR.rds",
#                     model_id,
#                     dataset_id,
#                     fold_id)
# 
# saveRDS(cond_LOO_SR,
#         filename)
# 
# ### Memory purge
# 
# rm(cond_LOO_pred,
#    cond_LOO_SR)

## Conditional - leave one in

### Load Data

filename <- sprintf("outputs/predictions/%s_%s_fold%s_condLOI.rds",
                    model_id,
                    dataset_id,
                    fold_id)

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

filename <- sprintf("outputs/species_richness/%s_%s_fold%s_condLOI_SR.rds",
                    model_id,
                    dataset_id,
                    fold_id)

saveRDS(cond_LOI_SR,
        filename)

### Memory purge

rm(cond_LOI_pred,
   cond_LOI_SR)

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

