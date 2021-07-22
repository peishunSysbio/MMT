################################################################
# To-do:  metabolomics analysis
# Author: Peishun Li, peishun@chalmers.se
################################################################
BARIA_wd <- "~/Documents/BARIA/baria"
outputDir <- "~/Documents/BARIA/baria/output"

####################################################
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
##
library(ggplot2)
library(DMwR)
library(reshape)
library(ggpubr)
library(ggsci)

##### Metabolomics data 
metabPeri201 <- read.table(paste(BARIA_wd,'/data/metabolomics/metabPeri201.txt',sep=""), header=T, sep="\t",stringsAsFactors=F)
identical(colnames(metabOriPeri201),colnames(metabPeri201))
subdesign201 <- read.table(paste(BARIA_wd,'/data/metabolomics/sub_design_Peripheral_blood201.txt',sep=""), header=T, sep="\t",stringsAsFactors=F, quote = "")
rownames(subdesign201)<-subdesign201$sampleId
# metabAnnoPeri201 <- read.table(paste(BARIA_wd,'/data/metabolomics/metabAnnoPeri201.txt',sep=""), header=T, sep="\t",stringsAsFactors=F, quote = "")
metabPeri201Norm <- as.data.frame(t(metabPeri201)/colSums(metabPeri201,na=T)*min(colSums(metabPeri201,na=T)))
identical(rownames(metabPeri201Norm),subdesign201$sampleId)

########### 1. two-way ANOVA with repeated measures using the aov() function in R to consider T2D status and time during MMTT (two factors)
identical(rownames(metabPeri201Norm),subdesign201$sampleId)
subdesign190 <- subset(subdesign201,clientId %in% as.numeric(rownames(groupInfo95)))
subdesign190$time <- ifelse(subdesign190$BLOOD.TYPE=='Peripheral blood', 'fasting', '2h')
subdesign190$time <- factor(subdesign190$time, levels = c('fasting', '2h'))
df <- groupInfo95[,c('t2d','groupFinal','Age', 'bmi', 'sex')]
df$clientId <- row.names(groupInfo95)
df <- merge(subdesign190, df, by='clientId')
rownames(df) <- df$sampleId
subdesign190 <- df[subdesign190$sampleId,]
subdesign190 <- rbind(subdesign190[1:95,][order(subdesign190$clientId[1:95]),],
subdesign190[96:190,][order(subdesign190$clientId[96:190]),])
#subdesign190$groupFinal <- subdesign190$groupFinal
#
sum(subdesign190$clientId[96:190]- as.numeric(rownames(groupInfo95)))
##
metabPeri190Norm <- metabPeri201Norm[subdesign190$sampleId,]
identical(rownames(metabPeri190Norm),subdesign190$sampleId)
metabPeri190LogNorm <- log(metabPeri190Norm)
metabAnno$metID
metabAnno$SUPER.PATHWAY[metabAnno$SUPER.PATHWAY==''] <-'Unknown'
metabAnno$SUB.PATHWAY[metabAnno$SUB.PATHWAY==''] <-'Unknown'
metabPeri190Anno <- metabAnno[colnames(metabPeri190LogNorm), ] 

df <- subdesign190[,c('groupFinal','clientId','time','Age', 'bmi', 'sex')]
df$clientId=factor(paste0('id',df$clientId))
table(df$sex)
#
aovTest <- data.frame()
for (i in 1:ncol(metabPeri190LogNorm)){
df$met <- metabPeri190LogNorm[,i]
# no adjust
#aov.out <-  summary(aov(met ~  groupFinal * time + Error(clientId/time), data=df))
# adjust for age
aov.out <-  summary(aov(met ~  Age + groupFinal * time + Error(clientId/time), data=df))
#
df1 <- cbind(as.data.frame(t(aov.out$`Error: clientId`[[1]]))[,c(1,2)], 
             as.data.frame(t(aov.out$`Error: clientId:time`[[1]]))[,c(1,2)])
#df1 <- cbind(as.data.frame(t(aov.out$`Error: clientId`[[1]])), 
 #            as.data.frame(t(aov.out$`Error: clientId:time`[[1]])))
df2 <- df1['Pr(>F)',c('groupFinal', 'time           ', 'groupFinal:time')]
#df2 <- df1['Pr(>F)',c('groupFinal', 'time           ', 'groupFinal:time')]
aovTest <- rbind(aovTest,df2)
}
rownames(aovTest) <- colnames(metabPeri190Norm)
colnames(aovTest) <- c('group_p', 'time_p', 'group_time_p')
sum(aovTest$group_p<0.05) 
sum(aovTest$time_p<0.05) 
sum(aovTest$group_time_p<0.05) 
## adjust p value by FDR
aovTest$group_pAdj <- p.adjust(aovTest$group_p, method="fdr")
aovTest$time_pAdj <- p.adjust(aovTest$time_p, method="fdr")
aovTest$group_time_pAdj <- p.adjust(aovTest$group_time_p, method="fdr")

#
aovTestSig <- cbind(aovTest, (aovTest < 0.05))
identical(metabPeri190Anno$metID, rownames(aovTestSig))
#
aovTestSig <- cbind(aovTestSig, metabPeri190Anno[,c('BIOCHEMICAL','SUPER.PATHWAY', 'SUB.PATHWAY')])
write.table(aovTestSig, file= paste(outputDir,'/MMTT_mets/aovTestSig_age_adj.txt',sep=""), sep="\t", quote=F,row.names = T)

#### dif metabolites among groups
sum(aovTest$group_pAdj<0.05) 
aovTestGroupDifmets <- rownames(aovTest)[aovTest$group_pAdj<0.05]
aovTestGroupDifmets145 <-aovTestGroupDifmets
intersect(aovTestGroupDifmets, aovTestGroupDifmets145)
#
aovTestGroupDifmets_adj <- rownames(aovTest)[aovTest$group_pAdj<0.05]
sum(!aovTestGroupDifmets %in% aovTestGroupDifmets_adj)
#
metabPeri95NormAnnoAovTestGroup <- cbind(metabAnno[aovTestGroupDifmets_adj, ],aovTestSig[aovTestGroupDifmets_adj,])

## Clasify differential metabolites into 3 types of response patterns
metabPeri95NormAnnoAovTestGroup$type <- 1
metabPeri95NormAnnoAovTestGroup$type <- ifelse(metabPeri95NormAnnoAovTestGroup[,22]=='FALSE'& metabPeri95NormAnnoAovTestGroup[,23]=='FALSE',1,
             ifelse(metabPeri95NormAnnoAovTestGroup[,22]=='TRUE'& metabPeri95NormAnnoAovTestGroup[,23]=='FALSE',2,3))
type3Mets <-  aovTestGroupDifmets_adj[metabPeri95NormAnnoAovTestGroup$type==3]
write.table(metabPeri95NormAnnoAovTestGroup, file= paste(outputDir,'/MMTT_mets/metabPeri95Anno_AovTestGroup145.txt',sep=""), sep="\t", quote=F,row.names = T)
table(metabPeri95NormAnnoAovTestGroup$type) 


#### metabolites changed over time
sum(aovTest$time_pAdj<0.05)  #439
aovTestTimeDifmets <- rownames(aovTest)[aovTest$time_pAdj<0.05]
metabPeri95NormAnnoAovTestTime <- metabAnno[aovTestTimeDifmets, ] 
metabPeri95LogNorMets439_timeDif <- metabPeri190LogNorm[, aovTestTimeDifmets]
table(metabPeri95NormAnnoAovTestTime$SUPER.PATHWAY)


########################### multiple pairwise comparisons
identical(rownames(metabPeri190Norm),subdesign190$sampleId)
df <- subdesign190[,c('groupFinal','clientId','time')]
df$group <- paste(subdesign190$groupFinal,subdesign190$time,sep='_')
cntrst <- c("preDiab_fasting-NGT_fasting", "T2D_fasting-NGT_fasting", "T2D_fasting-preDiab_fasting", 
            "preDiab_2h-NGT_2h", "T2D_2h-NGT_2h", "T2D_2h-preDiab_2h",
            "NGT_2h-NGT_fasting", "preDiab_2h-preDiab_fasting", "T2D_2h-T2D_fasting")
## stat info
source(paste0(BARIA_wd,"/src/diff_pairwise_compare.R", sep=""))
source(paste0(BARIA_wd,"/src/stat_group_outline.R", sep=""))
#
metabPeri190NormFC<- stat_group_outline(metabPeri190Norm, group= df$group,cntrst)
#
metabPeri190NormFC145mets <- metabPeri190NormFC[aovTestGroupDifmets_adj,]
row.names(metabPeri190NormFC145mets) <- metabPeri95NormAnnoAovTestGroup$BIOCHEMICAL


##### t tests for metabolomics data of 145 dif mets
metabPeri95LogNorMets145 <- metabPeri190LogNorm[, aovTestGroupDifmets_adj]
colnames(metabPeri95LogNorMets145) <- metabPeri95NormAnnoAovTestGroup$BIOCHEMICAL
identical(row.names(subdesign190), row.names(metabPeri95LogNorMets145))
##
metabPeri190Norm_ttest_group <- diff_pairwise_compare(metabPeri95LogNorMets145, df$group ,cntrst[1:6], difMethod='t.test', paired=F)

## diff mets among groups
metabPeri190Norm_ttest_group_sign <- data.frame(metabPeri190Norm_ttest_group[,c(2,4,6, 8,10,12)]<0.05) 
metabPeri190Norm_ttest_group_sign <- cbind(metabPeri190Norm_ttest_group_sign, metabPeri95NormAnnoAovTestGroup[,c('SUPER.PATHWAY', 'SUB.PATHWAY', 'BIOCHEMICAL')])
cntrst1 <- cntrst[1:6]
colnames(metabPeri190Norm_ttest_group_sign)[1:6] <- cntrst1
df<-data.frame(group=factor(cntrst1, level=cntrst1), mets=rep('t.test',6))
df$num<-colSums(metabPeri190Norm_ttest_group_sign[,1:6])
df$group1 <- rep(c('NGT-preDiab', 'NGT-T2D', 'preDiab-T2D'),2)
df$group2 <- factor(rep(c('fasting', '2h'),each=3), levels = c('fasting','2h'))

##
intersect(intersect(rownames(metabPeri190Norm_ttest_group_sign)[metabPeri190Norm_ttest_group_sign$`preDiab_2h-NGT_2h`],
rownames(metabPeri190Norm_ttest_group_sign)[metabPeri190Norm_ttest_group_sign$`T2D_2h-preDiab_2h`]),
rownames(metabPeri190Norm_ttest_group_sign)[metabPeri190Norm_ttest_group_sign$`T2D_2h-NGT_2h`])


###### calculae mets mean abudance
sampFile <- data.frame(groupFinal=subdesign190$groupFinal, time=subdesign190$time,row.names = row.names(subdesign190))
sampFile$group <- paste0(subdesign190$groupFinal, subdesign190$time)
#sampFile <- sampFile[,-c(1,2)]

#mat_t = merge(sampFile[3], metabPeri190Norm, by="row.names")
mat_t = merge(sampFile[3], log(metabPeri190Norm), by="row.names")

mat_t = mat_t[,-1]
mat_mean = aggregate(mat_t[,-1], by=mat_t[1], FUN=mean) 
groupName <- mat_mean[,1]
mat_mean_final <- t(mat_mean[,-1])
# mat_mean_final = do.call(rbind, mat_mean)[-1,]
colnames(mat_mean_final) = groupName
mat_mean_final145 <- mat_mean_final[aovTestGroupDifmets_adj,]
identical(rownames(mat_mean_final145), metabPeri95NormAnnoAovTestGroup$metID)
rownames(mat_mean_final145) <- metabPeri95NormAnnoAovTestGroup$BIOCHEMICAL
mat_mean_final145 <- mat_mean_final145[,c("NGTfasting", "NGT2h", "preDiabfasting", "preDiab2h", "T2Dfasting", "T2D2h")]

#  mets_mean  <- mat_mean_final

df <-  cbind(mets_mean, mets_sd, mets_cv)[aovTestGroupDifmets145,]
identical(aovTestGroupDifmets145, metabPeri95NormAnnoAovTestGroup$metID)
rownames(df) <- metabPeri95NormAnnoAovTestGroup$BIOCHEMICAL
write.table(df, file= paste(outputDir,'/MMTT_mets/met145-mean-sd.txt',sep=""), sep="\t", quote=F,row.names = T)


# Plotting the heatmap of dif met abundances
anno_col <- data.frame(timepoint=rep(c('fasting', '2h'),3), group = rep(c('NGT', 'preDiab', 'T2D'),each=2))
row.names(anno_col) <- colnames(mat_mean_final145)
anno_row <- data.frame(Pathway=metabPeri95NormAnnoAovTestGroup$SUPER.PATHWAY,
Type = as.factor(metabPeri95NormAnnoAovTestGroup$type), 
row.names = metabPeri95NormAnnoAovTestGroup$BIOCHEMICAL)
table(metabPeri95NormAnnoAovTestGroup$type)

anno_row$Pathway <- factor(anno_row$Pathway, levels = (c(as.character(unique(anno_row$Pathway)[1:6]),'Peptide','Unknown')))
anno_row <- anno_row[order(anno_row$Pathway),]
anno_row <- rbind(anno_row[anno_row$Type==1,],anno_row[anno_row$Type==2,],anno_row[anno_row$Type==3,])
mat_mean_final145sort <- mat_mean_final145[rownames(anno_row),]

##
bk <- c(seq(-2,-0.01,by=0.005),seq(0,2,by=0.005))
mycolor <- c(colorRampPalette(colors = c("#2c7bb6","white"))(length(bk)/2),colorRampPalette(colors = c("white","#d7191c"))(length(bk)/2))
anno_col$timepoint <- factor(anno_col$timepoint, levels = c('fasting','2h'))
ann_colors = list( timepoint = c('fasting'="#f0f0f0", '2h'="#bdbdbd"), 
                   group = c(NGT = "#2c7bb6", preDiab = "#fdae61", T2D='#d7191c'),
                   Pathway=c('Lipid'='#b3e2cd','Carbohydrate'='#fdcdac','Amino Acid'='#cbd5e8','Cofactors and Vitamins'='#f4cae4',
                             Nucleotide='#e6f5c9',Xenobiotics='#fff2ae',Peptide='#f1e2cc',Unknown='#737373'))
# 
p <- pheatmap(t(scale(t(mat_mean_final145sort))),
        cluster_rows=F, cluster_cols=F, main="",  fontsize = 7, fontsize_row=7, fontsize_col =7,
         annotation_col=anno_col,annotation_colors= ann_colors, annotation_row=anno_row1, 
         show_rownames=F, show_colnames =F, border_color=NA,cellwidth = 16, cellheight = 1.8,
         color = mycolor,legend_breaks=seq(-1.5,1.5, by = 0.5), breaks = bk,  annotation_legend=T,legend=T,
         gaps_col=c(2,4),gaps_row=c(sum(anno_row$Type==1),sum(anno_row$Type==1)+sum(anno_row$Type==2)))
# treeheight_row=10, treeheight_col=10, cellheight =10
ggsave(paste(outputDir,'/MMTT_mets/145difmets_heatmap1.pdf',sep=""), p, width = 6, height = 10)


# Calculating distances among NGT, Pre-D and T2D groups
sampleDist <- dist(t(mat_mean_final))
# sampleDist <- dist(scale(t(mat_mean_final)))
sampleDistMtx <- as.matrix(sampleDist)
# Plotting the distances: Group distances
sampleDistMtx <- sampleDistMtx[row.names(anno_col),row.names(anno_col)]
# Modify ordering of the clusters using clustering callback option
callback = function(hc, mat){
  sv = svd(t(mat))$v[,1]
  dend = reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}
anno_col$timepoint <- factor(anno_col$timepoint, levels = c('fasting','2h'))
ann_colors = list( timepoint = c('fasting'="#f0f0f0", '2h'="#bdbdbd"), 
                   group = c(NGT = "#2c7bb6", preDiab = "#fdae61", T2D='#d7191c'))
# 
p <- pheatmap(sampleDistMtx, cluster_rows=T, cluster_cols=T, main="",  fontsize = 7, fontsize_row=7, fontsize_col =7,
         annotation_col=anno_col, annotation_colors= ann_colors, show_rownames=F, show_colnames =F, treeheight_row=5, treeheight_col=5,
         clustering_callback = callback, cutree_rows=2, cutree_cols=2,annotation_legend=T,legend=T,
         cellwidth = 16, cellheight = 16, color=colorRampPalette(c("white","#ffeda0", "#756bb1"))(50)) # height = 10, width = 4
ggsave(paste(outputDir,'/MMTT_mets/group_distance_heatmap1.pdf',sep=""), p, width = 4, height = 4)



# Calculating distances among samples
sampleDist <- dist(log(metabPeri190Norm))
# sampleDist <- dist(scale(t(mat_mean_final)))
sampleDistMtx <- as.matrix(sampleDist)
identical(row.names(subdesign190), row.names(sampleDistMtx))
#
anno_col1 <- data.frame(timepoint=subdesign190$time, group = subdesign190$groupFinal)
row.names(anno_col1) <- rownames(subdesign190)

#anno_col1$timepoint <- factor(subdesign190$time, levels = c('fasting','2h'))
ann_colors1 = list( timepoint = c('fasting'="#f0f0f0", '2h'="#bdbdbd"), 
                   group = c(NGT = "#2c7bb6", preDiab = "#fdae61", T2D='#d7191c'))
# show_col( "#2c7bb6")
p <- pheatmap(sampleDistMtx, cluster_rows=T, cluster_cols=T, main="",  fontsize = 7, fontsize_row=7, fontsize_col =7,
              annotation_col=anno_col1, annotation_colors= ann_colors1, show_rownames=F, show_colnames =F, treeheight_row=5, treeheight_col=5,
               annotation_legend=T,legend=T, color=colorRampPalette(c("white","#ffeda0", "#756bb1"))(50),
              clustering_method='complete') # height = 10, width = 4, cellwidth = 10, cellheight = 10 
ggsave(paste(outputDir,'/MMTT_mets/samples_distance_heatmap1.pdf',sep=""), p, width = 4, height = 4)


