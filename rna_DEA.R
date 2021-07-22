################################################################
# To-do: differential expression analysis 
# Author: Peishun Li, peishun@chalmers.se
################################################################
BARIA_wd <- "~/Documents/BARIA/baria"
outputDir <- "~/Documents/BARIA/baria/output"

############
library(reshape)
library(ggpubr)
library(DESeq2)
################### Differential expression analysis between NGT, Pre-D and T2D groups
################################## Liver
# Creating the DESeqDataSet
#dds <- DESeqDataSetFromMatrix(countData=countLiver, colData=covars, design =~group) 
## Model adjusted by covariates/ confounding variable: age, bmi, sex
dds <- DESeqDataSetFromMatrix(countData=countLiver, colData=covars, design =~Age + bmi + sex + group)
################################## Jejunum
#dds <- DESeqDataSetFromMatrix(countData=countJejunum, colData=covars1, design =~group)
dds <- DESeqDataSetFromMatrix(countData=countJejunum, colData=covars1, design =~Age + bmi + sex + group)

#######
dds <-  estimateSizeFactors(dds)
sizeFactors(dds)
## deseq2 filter the low counts
## Only genes with the sum of counts across all smmples ≥10 
## and existed in at least five samples were considered in the analysis.
keep <-  (rowSums( counts(dds, normalized=TRUE) != 0 ) >= 0.05*ncol(dds)) & (rowSums(counts(dds, normalized=TRUE)) >= 10)
#
dds <- dds[keep,]
########################### DE test
dds <- DESeq(dds) 
dds
colnames(mcols(dds))
colnames(colData(dds)) 
resultsNames(dds)

##################### Extracting the DEA results
## T2D vs NGT
resT2D_NGT <- as.data.frame(results(dds, contrast=c("group","T2D","NGT")))
summary(resT2D_NGT)
colnames(resT2D_NGT) <- paste0(colnames(resT2D_NGT),'_T2D_NGT')
sum(resT2D_NGT$padj<0.05, na.rm=TRUE) 

## T2D vs preDiab
resT2D_PreD <- as.data.frame(results(dds, contrast=c("group","T2D" ,"preDiab")))
colnames(resT2D_PreD) <- paste0(colnames(resT2D_PreD),'_T2D_PreD')
sum(resT2D_PreD$padj<0.05, na.rm=TRUE)

## preDiab vs NGT
resPreD_NGT <- as.data.frame(results(dds, contrast=c("group","preDiab","NGT")))
colnames(resPreD_NGT) <- paste0(colnames(resPreD_NGT),'_PreD_NGT')
sum(resPreD_NGT$padj<0.05, na.rm=TRUE) 

## Liver-adjusted by covariates
resLiverAdj <- cbind(resPreD_NGT, resT2D_NGT, resT2D_PreD) 
colSums((resLiverAdj[c(6,12,18)]<0.05), na.rm = T)

## Jejunum-adjusted by covariates
resJejunumAdj <- cbind(resPreD_NGT, resT2D_NGT, resT2D_PreD)   
colSums((resJejunumAdj[c(6,12,18)]<0.05), na.rm = T)

###
df <- resLiverAdj
df <- resJejunumAdj

## all dif genes
sigGenesLiver <- rownames(df)[rowSums(df[,c(6,12,18)]<0.05, na.rm = T)>0]
sigGenesJejunum <- rownames(df)[rowSums(df[,c(6,12,18)]<0.05, na.rm = T)>0]
sigGenesVFat <- rownames(df)[rowSums(df[,c(6,12,18)]<0.05, na.rm = T)>0]
sigGenesSFat <- rownames(df)[rowSums(df[,c(6,12,18)]<0.05, na.rm = T)>0]

# overlap of dif tissues
# adj p<0.05
p <- venn.diagram(x= list( Liver = sigGenesLiver, Jejunum = sigGenesJejunum, VFat = sigGenesVFat, SFat =sigGenesSFat),
                  filename = NULL, height = 600, width = 600,resolution =300, imagetype="tiff",
                  fill=c('cornflowerblue', 'green', 'darkorchid1','orange'),
                  col="black", lwd = 1, alpha = 0.50, cex=0.45, cat.dist=0.04, cat.cex=0.45)

pdf(file=paste(outputDir, "/RNAseq/diffGeneVennAll.pdf", sep=""),width=1.6,height = 1.6)
##ggsave(file=paste(outputDir, "/RNAseq/diffGeneVennAll2.pdf", sep=""),width=4,height = 4, units ="cm")
grid.draw(p)
dev.off()


