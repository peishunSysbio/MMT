################################################################
# To-do: KEGG Gene Set Analysis of microbial KOs
# Author: Peishun Li, peishun@chalmers.se
################################################################
## get annotation of all KOs
koAnn <- read.table('/Users/peishun/Downloads/kegg/KOAnn1.txt',header=F, sep="\t",stringsAsFactors=F, quote = "")
colnames(koAnn) <- c('KO', 'Symbol', 'Name')
uniNames <- unique(koAnn$KO)
koAnn <- koAnn[match(uniNames,koAnn$KO),]
rownames(koAnn) <- koAnn$KO
rownames(resKOAdj)[!rownames(resKOAdj) %in% koAnn$KO]

############### KO-adjusted by covariates
resKOAnn <-resKOAdj[, c(2,5,6,8,11,12,14,17,18)]
resKOAnn <- merge(resKOAnn, koAnn, by="row.names", all.x=T)
sum(is.na(resKOAnn$KO))  
rownames(resKOAnn) <- resKOAnn$Row.names
resKOAnn <- resKOAnn[-1]

############################# use piano to do Gene Set Analysis of KEGG
## get KO-pathway maps
ko2pathway <- read.table('/Users/peishun/Downloads/kegg/ko2pathway1.txt',header=F, sep="\t",stringsAsFactors=F, quote = "")
colnames(ko2pathway) <- c('KEGG', 'KO')

## get pathway annotation
keggAnn <- read.table('/Users/peishun/Downloads/kegg/pathAnn1.txt',header=F, sep="\t",stringsAsFactors=F, quote = "")
colnames(keggAnn) <- c("KEGG", "Description")

##
KEGG_pathway <- merge(ko2pathway, keggAnn, by = "KEGG", all.x = TRUE)
KEGG_pathway <- KEGG_pathway[-1]
KEGG_pathway <- na.omit(KEGG_pathway) 
head(KEGG_pathway)
dim(unique(KEGG_pathway)) 

###### module
## get KO-module maps
ko2module <- read.table('/Users/peishun/Downloads/kegg/ko2module1.txt',header=F, sep="\t",stringsAsFactors=F, quote = "")
colnames(ko2module) <- c('module', 'KO')  ## 3158
temp <- as.data.frame(table(paste0(ko2module$module,ko2module$KO)))  ## 3088
ko2module <- unique(ko2module)

## get module annotation
moduleAnn <- read.table('/Users/peishun/Downloads/kegg/moduleAnn1.txt',header=F, sep="\t",stringsAsFactors=F, quote = "")
colnames(moduleAnn) <- c("module", "Description")
##
KEGG_module <- merge(ko2module, moduleAnn, by = "module", all.x = TRUE) 
KEGG_module <- KEGG_module[-1]
KEGG_module <- na.omit(KEGG_module) 
dim(unique(KEGG_module))  
head(KEGG_module)

## GSA
library(piano)
myGsc <- loadGSC(KEGG_pathway)  
##
myGsc <- loadGSC(KEGG_module)  

### load data for analysis
### choose p value and fold change
###  PreD vs NGT
allGenelist_PreD_NGT <- resKOAnn[,c(1:3,10:12)]
df <- allGenelist_PreD_NGT[!is.na(allGenelist_PreD_NGT$KO),]  ## 5641

### T2D vs NGT
allGenelist_T2D_NGT <- resKOAnn[,c(4:6,10:12)]
df <- allGenelist_T2D_NGT[!is.na(allGenelist_T2D_NGT$KO),]

###  T2D vs PreD
allGenelist_T2D_PreD <- resKOAnn[,c(7:9,10:12)]
df <- allGenelist_T2D_PreD[!is.na(allGenelist_T2D_PreD$KO),] 

###
##rownames(df) <- df$KO
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
resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7)
names(resList) <- c( "reporter1","reporter2","mean","median","sum","stouffer", "tailStrength")

##############################
### Save results as excel file, this table is important to see
# GSAsummaryTable(gsaRes1, save=TRUE, file=paste0(outputDir,'/microbiota/gsaResKegg_T2D_NGT_Liver.xls')) ##save file
gsaResTab <- GSAsummaryTable(gsaRes1, save=FALSE)
# To get info on a specific gene set
gs <- geneSetSummary(gsaRes1, "Phenylalanine metabolism")
gs$stats
### draw network plot
networkPlot(gsaRes1,class="non", adjusted = F, significance = 0.05, geneSets = NULL, overlap = 1, lay = 5,
            label = "names", cexLabel = 0.7, ncharLabel = 100, cexLegend = 1,nodeSize = c(5, 20),
            edgeWidth = c(1, 5), edgeColor = NULL, scoreColors = NULL, main= "non-directional")

#
glist <- resKOAnn[resKOAnn$KO %in% names(gs$geneLevelStats),]
glist$pvalue_PreD_NGT
glist[glist$pvalue_PreD_NGT<0.01,]

# Volcanol plot of fold change vs abundance plot
library(ggrepel)
# geom_text_repel, geom_label_repel
# T2D_NGT
glist1 <- glist[,c(4:6,10:12)]
glist1$significant <- 'no'
glist1$significant[glist1$pvalue_T2D_NGT < 0.01 & glist1$log2FoldChange_T2D_NGT < (-1)] <- 'down'
glist1$significant[glist1$pvalue_T2D_NGT < 0.01 & glist1$log2FoldChange_T2D_NGT > 1] <- 'up'

table(glist1$significant)
glist1$label <- glist1$KO
glist1$label[glist1$significant == 'no'] <- ""
summary(glist1$log2FoldChange_T2D_PreD)
#
p <- ggplot(glist1, aes(x=log2FoldChange_T2D_PreD, y=(-log10(pvalue_T2D_PreD)), color=significant, label = label)) + geom_point(size=0.5)  + 
  labs(x="log2(Fold Change)",y="-log10(P)")+main_theme+ geom_vline(xintercept=(-1),linetype=4,colour='grey')+ xlim(-4,4)+
  geom_vline(xintercept=(1),linetype=4,colour='grey')+      # values =c("#2c7bb6", "grey", "#d7191c" )
  geom_hline(yintercept=(-log10(0.01)),linetype=4,colour='grey')+ scale_color_manual(values =c("grey", "#d7191c"))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1),
        legend.position='none',legend.key.size = unit(0.2, "cm")) +
  geom_text_repel(size=2, segment.size = 0.2,segment.alpha =0.4) 
# geom_text(aes(label=rownames(df)), size=1,vjust = -0.5, nudge_y = 1) 
ggsave(paste(outputDir,'/microbiota/GSA/gs_path_T2D_NGT_Phenylalanine metabolism.pdf',sep=""), p, width = 4.5, height = 4, useDingbats=FALSE, units = "cm")

# heatmap of consensus Scores
pdf(paste(outputDir,'/microbiota/GSA/consensusScoresKegg_PreD_NGT.pdf',sep=""), width = 7, height = 8, onefile=FALSE)
ch <- consensusHeatmap(resList,cutoff=5,method="median",adjusted=T,ncharLabel = 60, colorkey = T)
dev.off()

