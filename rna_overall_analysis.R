################################################################
# To-do: PCA of RNA data
# Author: Peishun Li, peishun@chalmers.se
################################################################
library(ggplot2)
library(reshape)
library(dplyr)  #%>%
main_theme <- theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=7),
                   legend.position="right",
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=7),
                   text=element_text(family="Arial", size=7))

########################################################################################
## import clinical data
impVars <- read.table(file=paste(outputDir,'/clinic/impVars106.txt',sep=""), header=T, sep="\t",stringsAsFactors=T,row.names=1)
impVars$groupFinal
# enVars_b <- enVars
# groupInfo106 <- read.table(paste0(outputDir,'/groupInfo106.txt'), header=T, sep="\t",stringsAsFactors=F)
enVars_b <- enVars[,apply(enVars,2,function(x) sum(is.na(x)|(x=='')))==0]
covars <- enVars_b
colnames(covars)[5] <- 'group'
##
identical(rownames(impVars),rnaMetaLiver$patientID)
table(impVars$ppi, rnaMetaLiver$Ppi_Status)

# step 1: import RNAseq data
## Liver 
rnaCountLiver <- read.table(paste0(BARIA_wd,'/data/RNAseq/Counts_Liver.txt'), header=T, sep="\t",stringsAsFactors=F, row.names = 1)
rnaMetaLiver <- read.table(paste0(BARIA_wd,'/data/RNAseq/coldataLiver.txt'), header=T, sep="\t",stringsAsFactors=F)
identical(colnames(rnaCountLiver),rnaMetaLiver$X)
rnaMetaLiver$patientID<- gsub('Liver','' ,rnaMetaLiver$X)
# rownames(rnaMetaLiver) <- rnaMetaLiver$patientID
##
sum(rowSums(rnaCountLiver)==0) 
countLiver <- rnaCountLiver[rowSums(rnaCountLiver)!=0, ]
colnames(countLiver) <- rnaMetaLiver$patientID
identical(colnames(countLiver),rownames(covars))
summary(countLiver[,1])
rnaMetaLiver$group <- impVars$groupFinal

## Jejunum
rnaCountJejunum <- read.table(paste0(BARIA_wd,'/data/RNAseq/Counts_Jejunum.txt'), header=T, sep="\t",stringsAsFactors=F, row.names = 1)
rnaMetaJejunum <- read.table(paste0(BARIA_wd,'/data/RNAseq/coldataJejunum.txt'), header=T, sep="\t",stringsAsFactors=F)
identical(colnames(rnaCountJejunum),rnaMetaJejunum$X)
rnaMetaJejunum$patientID<- gsub('Jejunum','' ,rnaMetaJejunum$X)
rownames(rnaMetaJejunum) <- rnaMetaJejunum$patientID
## 
countJejunum <- rnaCountJejunum[rowSums(rnaCountJejunum)!=0, ]
identical(colnames(countJejunum),rnaMetaJejunum$X)
colnames(countJejunum) <- rnaMetaJejunum$patientID
summary(countJejunum[,1])
covars1 <- covars[rownames(rnaMetaJejunum),]
identical(colnames(countJejunum),rownames(covars1))

## Visceral fat 
rnaCountVFat <- read.table(paste0(BARIA_wd,'/data/RNAseq/countsVFat.txt'), header=T, sep="\t",stringsAsFactors=F, row.names = 1)
rnaMetaVFat <- read.table(paste0(BARIA_wd,'/data/RNAseq/coldataVFat.txt'), header=T, sep="\t",stringsAsFactors=F)
identical(colnames(rnaCountVFat),rnaMetaVFat$X)
rnaMetaVFat$patientID<- gsub('vFat','' ,rnaMetaVFat$X)
rownames(rnaMetaVFat) <- rnaMetaVFat$patientID
## 
countVFat <- rnaCountVFat[rowSums(rnaCountVFat)!=0, ]
identical(colnames(countVFat),rnaMetaVFat$X)
colnames(countVFat) <- rnaMetaVFat$patientID
summary(countVFat[,1])
covars1 <- covars[rownames(rnaMetaVFat),]

## Subcutaneous fat 
rnaCountSFat <- read.table(paste0(BARIA_wd,'/data/RNAseq/countsSFat.txt'), header=T, sep="\t",stringsAsFactors=F, row.names = 1)
rnaMetaSFat <- read.table(paste0(BARIA_wd,'/data/RNAseq/coldataSFat.txt'), header=T, sep="\t",stringsAsFactors=F)
identical(colnames(rnaCountSFat),rnaMetaSFat$X)
rnaMetaSFat$patientID<- gsub('sFat','' ,rnaMetaSFat$X)
rownames(rnaMetaSFat) <- rnaMetaSFat$patientID
## 
countSFat <- rnaCountSFat[rowSums(rnaCountSFat)!=0, ]
identical(colnames(rnaCountSFat),rnaMetaSFat$X)
colnames(countSFat) <- rnaMetaSFat$patientID
summary(countSFat[,1])
covars1 <- covars[rownames(rnaMetaSFat),]

#####################################
countAll <- cbind(rnaCountLiver, rnaCountJejunum, rnaCountVFat, rnaCountSFat)  ## 60675
identical(colnames(rnaMetaLiver[,-9]), colnames(rnaMetaVFat))
metaInfAll <- rbind(rnaMetaLiver[,-9],rnaMetaJejunum[,-9], rnaMetaVFat, rnaMetaSFat)
metaInfAll$type <- gsub('\\d*','',metaInfAll$X)
rownames(metaInfAll) <- metaInfAll$X

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=countAll, colData=metaInfAll, design =~group) 
dds <- DESeq(dds) 
############  visualization or clustering 
########## perform custom transformation - log of normalized counts, used in present study
countAllNorm <-  as.data.frame(t(counts(dds, normalized=TRUE)))
countAllNormLog <- as.data.frame(log10(countAllNorm + 1))
summary(as.numeric(countAllNormLog[1,]))

##############
## input data
table(metaInfAll$type)
# PCA principal component analysis 
library(ggbiplot)
mat<-countAllNormLog
matFilter <- mat[,-which(apply(mat, 2, var)==0)]
solPCA <- prcomp(scale(matFilter)) 

## ordination/ordinate of samples
ordi_site  <- data.frame(x = solPCA$x[,1], y = solPCA$x[,2], type=metaInfAll$type, group=metaInfAll$group)
PCAImp <- summary(solPCA)$importance[,c(1,2,3)]
# plot PCo 1 and 2
p <- ggplot(ordi_site, aes(x=x, y=y, color=group, shape=type)) + geom_point(size=0.5)+ main_theme  +
  labs(x=paste("PC1 (", format(100 *PCAImp[2,1], digits=4), "%)", sep=""),
       y=paste("PC2 (", format(100 *PCAImp[2,2], digits=4), "%)", sep=""),
       title="") +scale_colour_manual(values=c("#2c7bb6", "#fdae61", "#d7191c"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position='right',legend.key.size = unit(0.2, "cm")) 
ggsave(paste(outputDir,'/RNAseq/rnaAllPCA.pdf',sep=""), p, width = 4, height = 4, useDingbats=FALSE, units = "cm")




