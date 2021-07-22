##############################
# To-do:  mantel test of multi-omics
# Author: Boyang JI, Peishun Li, peishun@chalmers.se
##############################
# data import, tidy
library(tidyverse)
library(stringr)
library(ade4)
library(vegan)
library(usedist)

#-----------------------------
# load data
#-----------------------------

# datset name
all_list_data <- c("MGX-Taxon", "MGX-KOs", "HTX-Liver", "HTX-Jejunum", "HTX-Fat-Visceral", "HTX-Fat-SQ", "HMB_T0", "HMB_T2", "HMB_Ratio")

num_data_studies <- length(all_list_data)

# dataset directory
all_data_directory <- "./data/"
all_data_file <- c("MGX_species.txt", "MGX_ko.txt", "HTX_liver.txt", "HTX_jejunum.txt", "HTX_vfat.txt", "HTX_sfat.txt", "HMB_0h_s106.txt", "HMB_2h_s95.txt", "HMB_ratio.txt")

all_data_fullpath <- paste(all_data_directory, all_data_file, sep = "")
all_data_fullpath

# most of dataset included 106 samples, while the metabolomics data at 2h only have 95 samples
# so we will first get the number of samples, and corresponding SampleID of patients.
num_samples <- c(106, 106, 106, 105, 104, 105, 106, 95, 95)

# get the two type of samples
# sample_all_106 <- read_tsv("./data/HMB_0h_s106.txt") %>%
#   select(SampleID) %>%
#   arrange() %>%
#   unlist()
# sample_all_95 <- read_tsv("./data/HMB_2h_s95.txt") %>%
#   select(SampleID) %>%
#   arrange() %>%
#   unlist()
#
# index_95_in_106 <- match(sample_all_95, sample_all_106)

######################################################################
######################################################################
# get all the data
######################################################################

list_all_dataset <- list()

for (i in 1:num_data_studies) {
  i_file <- all_data_fullpath[i]

  i_data <- read_tsv(i_file)

  # no any pre-processing, because peishun had normalized them
  list_all_dataset[[i]] <- i_data
}

# save the Rdata file into the disk to save the loading time
save(list_all_dataset, file = "./data_output/data_list_7.RData")
# load("./data_output/data_list_7.RData")

######################################################################
######################################################################
# calculate the bray-curtis dist for each study
######################################################################
#
# using all sample 106/95
list_all_bray_dist <- list()
list_SampleID <- list()
list_samples_size <- vector()

for (i in 1:num_data_studies) {
  i_data <- list_all_dataset[[i]]

  i_sample_size <- nrow(i_data)
  list_samples_size[i] <- i_sample_size

  # if (i_sample_size == 106) {
  #   # i_test <- i_data %>% arrange(match(SampleID, sample_all_106[c(6, 1:5, 7)]))
  #   # i_test <- i_data %>% arrange(match(SampleID, sample_all_95))
  #   i_data_new <- i_data %>% arrange(match(SampleID, sample_all_106))
  # }
  # else if (i_sample_size == 95) {
  #   i_data_new <- i_data %>% arrange(match(SampleID, sample_all_95))
  # }

  i_data_new <- i_data %>% arrange(SampleID)
  i_SampleID <- i_data_new %>%
    select(SampleID) %>%
    unlist()
  list_SampleID[[i]] <- i_SampleID

  i_data_raw <- i_data_new[, -1]

  i_bray <- vegdist(i_data_raw, method = "bray")
  list_all_bray_dist[[i]] <- i_bray
}

######################################################################
######################################################################
# calculate the mantel test
######################################################################

list_all_mantel_test <- list()
combinations <- 0

for (i in 1:(num_data_studies - 1)) {
  for (j in i:num_data_studies) {
    # i_study <- all_list_data[i]
    # j_study <- all_list_data[j]

    i_dist <- list_all_bray_dist[[i]]
    j_dist <- list_all_bray_dist[[j]]

    i_sample_size <- list_samples_size[i]
    j_sample_size <- list_samples_size[j]

    i_sampleID <- list_SampleID[[i]]
    j_sampleID <- list_SampleID[[j]]

    common_sampleID <- intersect(i_sampleID, j_sampleID)
    i_common_index <- match(common_sampleID, i_sampleID)
    j_common_index <- match(common_sampleID, j_sampleID)

    i_dist_new <- dist_subset(i_dist, i_common_index)
    j_dist_new <- dist_subset(j_dist, j_common_index)

    ij_mantel <- mantel.rtest(i_dist_new, j_dist_new, nrepet = 9999)

    combinations <- combinations + 1
    list_all_mantel_test[[combinations]] <- ij_mantel
  }
}

######################################################################
######################################################################
# write the results to the table
######################################################################

mantel_test_result <- tibble(
  Study1 = character(),
  Study2 = character(),
  Obs = numeric(),
  Pvalue = numeric(),
  StdObs = numeric(),
  Expection = numeric(),
  Variance = numeric()
)

# set a count number
out_row_count <- 0

for (i in 1:(num_data_studies - 1)) {
  for (j in i:num_data_studies) {
    i_study <- all_list_data[i]
    j_study <- all_list_data[j]
    
    out_row_count <- out_row_count + 1
    
    if (i == j) {
      next
    }
    else {
      i_mantel_data <- list_all_mantel_test[[out_row_count]]
      i_obs <- i_mantel_data$obs
      i_P <- i_mantel_data$pvalue
      i_std_obs <- i_mantel_data$expvar[1]
      i_exp <- i_mantel_data$expvar[2]
      i_var <- i_mantel_data$expvar[3]
    }
    
    # mantel_test_result <- rbind(mantel_test_result, c(i_study, j_study, i_obs, i_P, i_std_obs, i_exp, i_var ))
    mantel_test_result <- add_row(mantel_test_result, Study1=i_study, Study2=j_study, Obs=i_obs, Pvalue=i_P, StdObs=i_std_obs, Expection=i_exp, Variance=i_var )
  }
}

mantel_test_result
write_tsv(mantel_test_result, "./data_output/data_list_mantel_test_inter_individual.txt")

