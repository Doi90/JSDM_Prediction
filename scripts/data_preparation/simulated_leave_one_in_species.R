##########################################################
##########################################################
###                                                    ###
###  DETERMINE LEFT IN SPECIES FOR SIMULATED DATASETS  ###
###                                                    ###
###   This script determines the species that are left ###
### in for the conditional leave one in prediction     ###
### method. Saves to file as an RDS                    ###
###                                                    ###
##########################################################
##########################################################

datasets <- c("sim1random",
              "sim2random",
              "sim3random",
              "sim4random",
              "sim5random",
              "sim6random",
              "sim7random",
              "sim8random",
              "sim9random",
              "sim10random",
              "sim1spatial",
              "sim2spatial",
              "sim3spatial",
              "sim4spatial",
              "sim5spatial",
              "sim6spatial",
              "sim7spatial",
              "sim8spatial",
              "sim9spatial",
              "sim10spatial")

spp_LOI_storage <- data.frame(low = character(length(datasets)),
                              med = character(length(datasets)),
                              high = character(length(datasets)),
                              stringsAsFactors = FALSE)

rownames(spp_LOI_storage) <- datasets


for(i in datasets){
  
  tmp <- readRDS(sprintf("data/%1$s/%1$s_communitySimul.rds",
                         i))
  
  spp_richness <- colSums(tmp$data$Y)
  
  low_spp <- which(names(spp_richness) == names(spp_richness)[order(spp_richness)[1]])
  
  med_spp <- which(names(spp_richness) == names(spp_richness)[order(spp_richness)[5]])
  
  high_spp <- which(names(spp_richness) == names(spp_richness)[order(spp_richness)[10]])
  
  spp_LOI_storage[which(datasets == i), ] <- c(low_spp,
                                               med_spp,
                                               high_spp)
  
}

saveRDS(spp_LOI_storage,
        "data/simulated_datasets_species_left_in.rds")

for(i in datasets){
  
  tmp <- readRDS(sprintf("data/%1$s/%1$s_communitySimul.rds",
                         i))
  
  spp_richness <- colSums(tmp$data$Y)
  
  print(spp_richness)
  
}
