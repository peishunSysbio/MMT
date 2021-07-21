################################################################
# To-do: prepare data for predictive model using ridge regression
# Author: Peishun Li, peishun@chalmers.se
################################################################
library(ggplot2)
library(reshape)
library(ggpubr)
library(tidyverse)
BARIA_wd <- "~/Documents/BARIA/baria"
outputDir <- "~/Documents/BARIA/baria/output"

################################## load data
## Using the profiles of taxa, KOs, and metabolites
############################   meta data
data_dir <- '/Users/peishun/Documents/BARIA/baria/output/ML_regression'
metadata_file <- paste0(data_dir, '/omicsData20191115/impVars106.txt')
# metadata <- read_tsv(metadata_file)
metadata <- read.table(metadata_file, header=T, sep="\t",stringsAsFactors=F, quote = "",row.names=1)
table(metadata$group)
# change name
#metadata$groupFinal <- factor(metadata$groupFinal, labels = c("NGT", "Prediabetes", "T2D"))
##
group_info_useful <- metadata %>% select(Age, bmi, totalAUC_gluc)
## converting the categorical variables into dummy variables (0 and 1 numeric values) 
## recode your factor variables using dummy variables,  contrasts(factor(metadata$sex))
group_info_useful$sex1 <- model.matrix(~metadata$sex)[,2]
table(group_info_useful$sex1,metadata$sex)

##################################  microbiota species
data_file <- paste0(data_dir, '/omicsData20191115/metagenomics/speciePer.txt')
#microbiota_data <- read_tsv(microbiota_file)
microbiota_data <- read.table(data_file, header=T, sep="\t",stringsAsFactors=F, quote = "",row.names=1)
microbiota_species <- colnames(microbiota_data)
# 
col_names_new <- paste("s", 1:ncol(microbiota_data), sep = "_")
colnames(microbiota_data) <- col_names_new
head(microbiota_data)
species_info <- data.frame(name = col_names_new, species= microbiota_species)
# log transformation
summary(unlist(microbiota_data_log))
microbiota_data_log <- log10(microbiota_data + 1e-10)

###################################  mets data
metabAnno <- read.table(paste0(data_dir, '/omicsData20191115/metabolomics/metabAnno.txt'), header=T, sep="\t",stringsAsFactors=F, quote = "")
rownames(metabAnno) <- metabAnno$metID
## fasting
metabPeri106Norm <- read.table(paste0(data_dir, '/omicsData20191115/metabolomics/metabPeri106Norm.txt',sep=""), header=T, sep="\t",stringsAsFactors=F,row.names=1)
identical(rownames(metabPeri106Norm),rownames(metadata))

# log transformation
summary(unlist(metabPeri106Norm_log))
metabPeri106Norm_log <- log10(metabPeri106Norm)
metabPeri106Norm_log <- metabPeri106Norm_log[, !colnames(metabPeri106Norm_log)=='met541']

## 2h
metabPeri2h95Norm <- read.table(paste0(data_dir, '/omicsData20191115/metabolomics/metabPeri2h95Norm.txt',sep=""), header=T, sep="\t",stringsAsFactors=F,row.names=1)
metadata95 <- metadata[rownames(metabPeri2h95Norm),]
identical(rownames(metabPeri2h95Norm),rownames(metadata95))
# log transformation
summary(unlist(metabPeri2h95Norm_log))
metabPeri2h95Norm_log <- log10(metabPeri2h95Norm)
metabPeri2h95Norm_log <- metabPeri2h95Norm_log[, !colnames(metabPeri2h95Norm_log)=='met541']

## ratio
metabPeri95NormRatio <- read.table(paste0(data_dir, '/omicsData20191115/metabolomics/metabPeri95NormRatio.txt',sep=""), header=T, sep="\t",stringsAsFactors=F,row.names=1)
identical(rownames(metabPeri95NormRatio),rownames(metadata95))
# log transformation
summary(unlist(metabPeri95NormRatio_log))
metabPeri95NormRatio_log <- log2(metabPeri95NormRatio)
metabPeri95NormRatio_log <- metabPeri95NormRatio_log[, !colnames(metabPeri95NormRatio_log)=='met541']


###################################  RNAseq data
rnaAnno <- read.table(paste0(data_dir, '/omicsData20191115/RNAseq/idAllGenes.txt'), header=T, sep="\t",stringsAsFactors=F, quote = "")

## Liver
countLiverNorm <- read.table(paste0(data_dir, '/omicsData20191115/RNAseq/countLiverNorm.txt',sep=""), header=T, sep="\t",stringsAsFactors=F,row.names=1)
identical(rownames(countLiverNorm),rownames(metadata))
# log transformation
summary(unlist(countLiverNorm_log))
countLiverNorm_log <- log10(countLiverNorm+1)

## Jejunum
countJejunumNorm <- read.table(paste0(data_dir, '/omicsData20191115/RNAseq/countJejunumNorm.txt',sep=""), header=T, sep="\t",stringsAsFactors=F,row.names=1)
identical(rownames(countJejunumNorm),rownames(metadata))
# log transformation
summary(unlist(countJejunumNorm_log))
countJejunumNorm_log <- log10(countJejunumNorm+1)

## VFat
countVFatNorm <- read.table(paste0(data_dir, '/omicsData20191115/RNAseq/countVFatNorm.txt',sep=""), header=T, sep="\t",stringsAsFactors=F,row.names=1)
identical(rownames(countVFatNorm),rownames(metadata))
# log transformation
summary(unlist(countVFatNorm_log))
countVFatNorm_log <- log10(countVFatNorm+1)

## SFat
countSFatNorm <- read.table(paste0(data_dir, '/omicsData20191115/RNAseq/countSFatNorm.txt',sep=""), header=T, sep="\t",stringsAsFactors=F,row.names=1)
identical(rownames(countSFatNorm),rownames(metadata))
# log transformation
summary(unlist(countSFatNorm_log))
countSFatNorm_log <- log10(countSFatNorm+1)

## group inf
group_info_useful95 <- group_info_useful[rownames(metabPeri2h95Norm), ]
#
group_info_Jejunum <- group_info_useful[rownames(countJejunumNorm_log), ]
#
group_info_VFat <- group_info_useful[rownames(countVFatNorm_log), ]
#
group_info_SFat <- group_info_useful[rownames(countSFatNorm_log), ]

##################################
## Microbiota
omicData <- microbiota_data_log


################################# Integrate data
## species + mets_0h
identical(row.names(metabPeri106Norm_log),row.names(microbiota_data_log))
metSpe <- cbind(metabPeri106Norm_log,microbiota_data_log)
omicData <- metSpe
ncol(microbiota_data_log)+ncol(metabPeri106Norm_log)
##
metSpeAnno <- rbind(setNames(metabAnno[,c(14,2),],names(species_info)),species_info)
colnames(metSpeAnno) <- c('ID', 'name')
rownames(metSpeAnno) <- metSpeAnno$ID


# all_data <- full_join(group_info_useful, omicData, by = "SampleID")
all_data <- cbind(group_info_useful, omicData)
#
all_data <- cbind(group_info_useful95, omicData)
#
all_data <- cbind(group_info_Jejunum, omicData)
#
head(all_data)
identical(rownames(metadata),row.names(all_data))
all_data <- all_data[metadata$groupFinal !='T2D',]  #84
#
X_data <- all_data %>% dplyr::select(-totalAUC_gluc) %>% as.matrix()
Y_data <- all_data %>% dplyr::select(totalAUC_gluc) %>% as.matrix()

