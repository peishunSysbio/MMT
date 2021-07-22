################################################################
# To-do:  Differential analysis of microbial species and KOs using DESeq2
# Author: Peishun Li, peishun@chalmers.se
################################################################
library(dplyr)
library(limma)
library(ggplot2)
library(reshape2)
library(ggpubr)
BARIA_wd <- "~/Documents/BARIA/baria"
outputDir <- "~/Documents/BARIA/baria/output"

###############################################################################################
#  Pairwise diff analysis among NGT, Per-D and T2D groups 
covars <- enVars_b
colnames(covars)[5] <- 'group'
# speces
counts <- specieCount
# or KOs
counts <- koCount

######## univariate generalized linear models(GLM)
#dds <- DESeqDataSetFromMatrix(countData=counts, colData=covars, design =~group)

######## Model adjusted by covariates: age, bmi, sex
dds <- DESeqDataSetFromMatrix(countData=counts, colData=covars, design =~Age + bmi + sex + group)

#######
dds <-  estimateSizeFactors(dds)
## deseq2 filter the low counts
## Only species with the sum of counts across all smmples ≥10 
## and existed in at least five samples were considered in the analysis.
keep <-  (rowSums( counts(dds, normalized=TRUE) != 0 ) >= 0.05*ncol(dds)) & (rowSums(counts(dds, normalized=TRUE)) >= 10)
##
dds <- dds[keep,]

## Performing estimation of dispersion parameter
dds.disp <- estimateDispersions(dds)
## A diagnostic plot which shows the mean of normalized counts (x axis)
## and dispersion estimate for each genes
plotDispEsts(dds.disp)
########################### DE test
dds <- DESeq(dds) 
dds
colnames(mcols(dds))

## preDiab vs NGT
resPreD_NGT <- as.data.frame(results(dds, contrast=c("group","preDiab","NGT")))
colnames(resPreD_NGT) <- paste0(colnames(resPreD_NGT),'_PreD_NGT')
# sum(resPreD_NGT$pvalue<0.05, na.rm=TRUE) 
sum(resPreD_NGT$padj<0.05, na.rm=TRUE) #26, covar-adj 13

## T2D vs NGT
resT2D_NGT <- as.data.frame(results(dds, contrast=c("group","T2D","NGT")))
summary(resT2D_NGT)
colnames(resT2D_NGT) <- paste0(colnames(resT2D_NGT),'_T2D_NGT')
sum(resT2D_NGT$pvalue<0.05, na.rm=TRUE)
sum(resT2D_NGT$padj<0.05, na.rm=TRUE) # 99, covar-adj 59
hist(resT2D_NGT$padj, breaks=20, col="grey", main="DESeq2 p-value distribution", xlab="DESeq2 P-value", ylab="Number of genes")

## T2D vs preDiab
resT2D_PreD <- as.data.frame(results(dds, contrast=c("group","T2D" ,"preDiab")))
colnames(resT2D_PreD) <- paste0(colnames(resT2D_PreD),'_T2D_PreD')
sum(resT2D_PreD$pvalue<0.05, na.rm=TRUE) 
sum(resT2D_PreD$padj<0.05, na.rm=TRUE)  #94, covar-adj 59
##
identical(rownames(resT2D_NGT), rownames(resPreD_NGT))

##  adjusted by covariates
resAdj <- cbind(resPreD_NGT, resT2D_NGT, resT2D_PreD)

## Draw HeatMap of speces fold changes
mulvar_pValue <- resAdj[impSpecies,c(6,12,18)]
colnames(mulvar_pValue) <- c("preDiab-NGT", "T2D-NGT",  "T2D-preDiab")
mulvar_pValue[is.na(mulvar_pValue)] <- 1
mulvar_fc <- resAdj[impSpecies,c(2,8,14)]
colnames(mulvar_fc) <- c("preDiab-NGT", "T2D-NGT",  "T2D-preDiab")
summary(unlist(mulvar_fc))
mulvar_fc1 <- mulvar_fc
mulvar_fc1[mulvar_fc1>5] <- 5
mulvar_fc1[mulvar_fc1<(-5)] <- (-5)
##
library(pheatmap)
# pheatmap(corMatrix, display_numbers = TRUE,cluster_rows  = FALSE,cluster_cols = FALSE) 
bk <- c(seq(-5,-0.01,by=0.05),seq(0,5,by=0.05))
mycolor <- c(colorRampPalette(colors = c("#2c7bb6","white"))(length(bk)/2),colorRampPalette(colors = c("white","#d7191c"))(length(bk)/2))
p <- pheatmap(t(mulvar_fc1), display_numbers = t(matrix(ifelse(mulvar_pValue < 0.01, "*", 
              ifelse(mulvar_pValue< 0.05, "+",  "")), nrow(mulvar_pValue))), scale='none',cellwidth=7,
              cluster_rows  = F,cluster_cols = T,fontsize_number=6,fontsize=6,silent=T,angle_col = c("90"), 
              color = mycolor, legend = F, legend_breaks=seq(-5,5,1), breaks = bk, border_color ='grey',treeheight_row=0,treeheight_col =0)  
ggsave(paste(outputDir,'/microbiota/imp_species_diff_fc_deseq1.pdf',sep=""), p, width = 6, height = 8, units="cm", useDingbats=FALSE)


