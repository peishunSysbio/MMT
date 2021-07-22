################################################################
# To-do: Gene Set Analysis of KEGG for different tissues
# Author: Peishun Li, peishun@chalmers.se
################################################################
BARIA_wd <- "~/Documents/BARIA/baria"
outputDir <- "~/Documents/BARIA/baria/output"

################################################## Gene Set Analysis
## get annotation of all genes
library(org.Hs.eg.db)
library(clusterProfiler)
keytypes(org.Hs.eg.db)
ids <- bitr(rownames(countAll), fromType="ENSEMBL", toType=c("ENTREZID", "SYMBOL","GENENAME"), OrgDb="org.Hs.eg.db")
uniNames <- unique(ids$ENSEMBL)
ids <- ids[match(uniNames,ids$ENSEMBL),]
rownames(ids) <- ids$ENSEMBL
ids <- ids[-1] ## 25478
idAllGenes <- data.frame(ENSEMBL=rownames(countAll))
rownames(idAllGenes) <- idAllGenes$ENSEMBL
idAllGenes <- merge(idAllGenes, ids, by="row.names", all.x=T)
sum(is.na(idAllGenes$ENTREZID))  
rownames(idAllGenes) <- idAllGenes$ENSEMBL
idAllGenes <- idAllGenes[-1]

############### Liver-adjusted by covariates
resLiverAdj <- read.table(paste0(outputDir,'/RNAseq/resLiverAdj.txt'), header=T, sep="\t",stringsAsFactors=F, row.names = 1)
resLiverAnn <-resLiverAdj[, c(2,5,6,8,11,12,14,17,18)]
resLiverAnn <- merge(resLiverAnn, idAllGenes, by="row.names", all.x=T)
sum(is.na(resLiverAnn$ENTREZID))  
rownames(resLiverAnn) <- resLiverAnn$Row.names
resLiverAnn <- resLiverAnn[-1]

############### Jejunum adjusted by covariates
resJejunumAdj <- read.table(paste0(outputDir,'/RNAseq/resJejunumAdj.txt'), header=T, sep="\t",stringsAsFactors=F, row.names = 1)
resJejunumAnn <-resJejunumAdj[, c(2,5,6,8,11,12,14,17,18)]
resJejunumAnn <- merge(resJejunumAnn, idAllGenes, by="row.names", all.x=T)
sum(is.na(resJejunumAnn$ENTREZID))  
rownames(resJejunumAnn) <- resJejunumAnn$Row.names
resJejunumAnn <- resJejunumAnn[-1]

############ VFat
resVFatAdj <- read.table(paste0(outputDir,'/RNAseq/resVFatAdj.txt'), header=T, sep="\t",stringsAsFactors=F, row.names = 1)
resVFatAnn <-resVFatAdj[, c(2,5,6,8,11,12,14,17,18)]
resVFatAnn <- merge(resVFatAnn, idAllGenes, by="row.names", all.x=T)
sum(is.na(resVFatAnn$ENTREZID))  
rownames(resVFatAnn) <- resVFatAnn$Row.names
resVFatAnn <- resVFatAnn[-1]

############ SFat
resSFatAdj <- read.table(paste0(outputDir,'/RNAseq/resSFatAdj.txt'), header=T, sep="\t",stringsAsFactors=F, row.names = 1)
resSFatAnn <-resSFatAdj[, c(2,5,6,8,11,12,14,17,18)]
resSFatAnn <- merge(resSFatAnn, idAllGenes, by="row.names", all.x=T)
sum(is.na(resSFatAnn$ENTREZID))  
rownames(resSFatAnn) <- resSFatAnn$Row.names
resSFatAnn <- resSFatAnn[-1]

#############################  use piano to do Gene Set Analysis
######### Gene Set Analysis of KEGG
# BiocManager::install("KEGG.db")
# library(KEGG.db)   #  This package should now be considered deprecated
library(piano)
library(clusterProfiler)
hsa_kegg <- download_KEGG('hsa')
# bitr_kegg()
res1 <- hsa_kegg$KEGGPATHID2EXTID
res2 <- hsa_kegg$KEGGPATHID2NAME
head(res1)
colnames(res1) <- c("KEGG", "ENTREZID")
colnames(res2) <- c("KEGG", "Description")
KEGG_pathway <- merge(res1, res2, by = "KEGG", all.x = TRUE)
KEGG_pathway <- KEGG_pathway[-1]
head(KEGG_pathway)
dim(unique(KEGG_pathway))  
myGsc <- loadGSC(KEGG_pathway)   

### load data for analysis
### choose significant and fold change
################## Liver
all(unique(KEGG_pathway$ENTREZID) %in% resLiverAnn$ENTREZID)
length(intersect(unique(KEGG_pathway$ENTREZID), resLiverAnn$ENTREZID))  

###  PreD vs NGT
allGenelist_PreD_NGT <- resLiverAnn[,c(1:3,10:12)]
df <- allGenelist_PreD_NGT[!is.na(allGenelist_PreD_NGT$ENTREZID),] 

### T2D vs NGT
allGenelist_T2D_NGT <- resLiverAnn[,c(4:6,10:12)]
df <- allGenelist_T2D_NGT[!is.na(allGenelist_T2D_NGT$ENTREZID),] 

###  T2D vs PreD
allGenelist_T2D_PreD <- resLiverAnn[,c(7:9,10:12)]
df <- allGenelist_T2D_PreD[!is.na(allGenelist_T2D_PreD$ENTREZID),] 

###
uniNames <- unique(df$ENTREZID)
df <- df[match(uniNames,df$ENTREZID),]
rownames(df) <- df$ENTREZID
head(df)
myPval <- df[2]
head(myPval)
myFC <- df[1]
head(myFC)

### run GSA
gsaRes1 <- runGSA(myPval, directions = myFC, gsc = myGsc, geneSetStat = "reporter", signifMethod="geneSampling",
                  gsSizeLim = c(10,300), nPerm=1000)
gsaRes2 <- runGSA(myPval, directions = myFC, gsc = myGsc, geneSetStat = "reporter", signifMethod="nullDist",
                  gsSizeLim = c(10,300), nPerm=1000)
gsaRes3 <- runGSA(myPval, directions = myFC, gsc = myGsc, geneSetStat = "mean", signifMethod="geneSampling",
                  gsSizeLim = c(10,300), nPerm=1000)
gsaRes4 <- runGSA(myPval, directions = myFC, gsc = myGsc, geneSetStat = "median", signifMethod="geneSampling",
                  gsSizeLim = c(10,300), nPerm=1000)
gsaRes5 <- runGSA(myPval, directions = myFC, gsc = myGsc, geneSetStat = "sum", signifMethod="geneSampling",
                  gsSizeLim = c(10,300), nPerm=1000)
gsaRes6 <- runGSA(myPval, directions = myFC, gsc = myGsc, geneSetStat = "stouffer", signifMethod="geneSampling",
                  gsSizeLim = c(10,300), nPerm=1000)
gsaRes7 <- runGSA(myPval, directions = myFC, gsc = myGsc, geneSetStat = "tailStrength", signifMethod="geneSampling",
                  gsSizeLim = c(10,300), nPerm=1000)
             gsSizeLim = c(10,300), nPerm=1000)
##############################
resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7)
names(resList) <- c( "reporter1","reporter2","mean","median","sum","stouffer", "tailStrength")

### Save results as excel file, this table is important to see
gsaResTab <- GSAsummaryTable(gsaRes1, save=FALSE)

# To get info on a specific gene set
gs <- geneSetSummary(gsaRes1, "Insulin secretion")
#
glist <- resLiverAnn[resLiverAnn$ENTREZID %in% names(gs$geneLevelStats),]
#
################
# Volcanol plot of fold change vs abundance plot
library(ggrepel)
# T2D_NGT
glist1 <- glist[,c(4:6,11:13)]
glist1$significant <- 'no'
glist1$significant[glist1$pvalue_T2D_NGT < 0.05 & glist1$log2FoldChange_T2D_NGT < (-0.2)] <- 'down'
glist1$significant[glist1$pvalue_T2D_NGT < 0.05 & glist1$log2FoldChange_T2D_NGT > 0.2] <- 'up'
table(glist1$significant)
glist1$label <- glist1$SYMBOL
glist1$label[glist1$significant == 'no'] <- ""
summary(glist1$log2FoldChange_T2D_NGT)
p <- ggplot(glist1, aes(x=log2FoldChange_T2D_NGT, y=(-log10(pvalue_T2D_NGT)), color=significant, label = label)) + geom_point(size=0.5)  + 
  labs(x="log2(Fold Change)",y="-log10(P)")+main_theme+ geom_vline(xintercept=(-0.2),linetype=4,colour='grey')+ xlim(-1,1)+
  geom_vline(xintercept=(0.2),linetype=4,colour='grey')+      # values =c("#2c7bb6", "grey", "#d7191c" )
  geom_hline(yintercept=(-log10(0.05)),linetype=4,colour='grey')+ scale_color_manual(values =c("#2c7bb6", "grey", "#d7191c"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position='none',legend.key.size = unit(0.2, "cm")) +
  geom_text_repel(size=2, segment.size = 0.2,segment.alpha =0.4) 
# geom_text(aes(label=rownames(df)), size=1,vjust = -0.5, nudge_y = 1) 

### draw network plot
networkPlot(gsaRes1,class="non", adjusted = T, significance = 0.05, geneSets = NULL, overlap = 1, lay = 5,
            label = "names", cexLabel = 0.7, ncharLabel = 100, cexLegend = 1,nodeSize = c(5, 20),
            edgeWidth = c(1, 5), edgeColor = NULL, scoreColors = NULL, main= "non-directional")

####  heatmap of consensus Scores
pdf(paste(outputDir,'/RNAseq/GSA/consensusScoresKegg_PreD_NGT_SFat_median_5.pdf',sep=""), width = 7, height = 8, onefile=FALSE)
ch <- consensusHeatmap(resList,cutoff=5,method="median",adjusted=T,ncharLabel = 60, colorkey = T)
dev.off()





