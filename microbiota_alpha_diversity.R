################################################################
# To-do:  alpha diversity analysis of microbiota
# Author: Peishun Li, peishun@chalmers.se
################################################################
library(vegan)
library(dplyr) 
library(ggplot2)

################################
# calculate the alpha diversity
################################
# relative abundance of species level
comm <- speciePer

# /Shannon/Simpson/invsimpson
## shannon Diversity Indices
shannonSpecie <- diversity(comm, index="shannon")
summary(shannonSpecie)
# simpson
simpSpecie <- diversity(comm, index="simpson")
summary(simpSpecie)
# invsimpson
invSpecie <- diversity(comm, index="invsimpson")
summary(invSpecie)

# species richness - observed number of species
richSpecie <- specnumber(comm)  ## rowSums(comm > 0) 

##############################
specieDiversityAlpha <- data.frame(shannon=shannonSpecie,simpson=simpSpecie,invsimpson=invSpecie, 
                                   species_richnes=richSpecie, group=microbDesign106$group, row.names=microbDesign106$clientId)
# write species level diversity to file
write.table(specieDiversityAlpha, file=paste0(outputDir,"/microbiota/alpha_species_rel_abu.txt"), sep="\t", quote=F)

############ statistical analysis of alpha diversity for the NGT, Pre-D and T2D groups
# anova
specieDiversity_stats <- aov(simpson ~ group, data = specieDiversityAlpha)   
summary(specieDiversity_stats)   

## kruskal test
kruskal.test(simpson ~ group, data = specieDiversityAlpha)



