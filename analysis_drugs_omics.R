################################################################
# To-do: analysis effects of drugs on multi-omics
# Author: Peishun Li, peishun@chalmers.se
################################################################
## import drug vars
drugVars <- read.table(file=paste(outputDir,'/clinic/drugVars.txt',sep=""), header=T, sep="\t",stringsAsFactors=T,row.names=1)

# species profile
comm <- speciePer
# metabolomics
comm <- log(metabPeri106Norm)
# Liver RNA
#countLiverNormLog <- as.data.frame(log10(countLiverNorm + 1))
comm <- countLiverNormLog

####################################### 
# check the potential confounders of the association between disease group and 
# gut microbiota composition maybe play an important role.
# import the potential confounders
enVars_b <- drugVars
enVars_b <- drugVars[rownames(comm),]
apply(drugVars,2,function(x) sum(is.na(x)|(x=='')))
identical(rownames(enVars_b),rownames(comm))
####################### PERMANOVA R2 and P value for univariate comparison
df <-data.frame()
for(j in 1:ncol(enVars_b)) {
  if(sum(is.na(enVars_b[,j]))>0){
  comm1 <- comm[-which(is.na(enVars_b[,j])),]}
  else{
    comm1 <- comm
  }
  set.seed(123)
  adonis_envs <- adonis(comm1 ~ ., data = na.omit(enVars_b[,j,drop=F]), permutations=999, method = "bray")
  df <- rbind(df, data.frame(envs=rownames(adonis_envs$aov.tab)[1], r2=adonis_envs$aov.tab$R2[1],
             pvals=adonis_envs$aov.tab$`Pr(>F)`[1]))
}
rownames(na.omit(enVars_b[,j,drop=F]))
df$p_adj <- p.adjust(df$pvals, method="fdr")
#
#adonis_envs_univariate <- df
# species
write.table(df, file= paste(outputDir,'/microbiota/adonis_drugs_univariate.txt',sep=""), sep="\t", quote=F,row.names = F)
# fasting metabolites 
write.table(df, file= paste(outputDir,'/MMTT_mets/adonis_drugs_univariate_mets0h.txt',sep=""), sep="\t", quote=F,row.names = F)
# Liver RNA
write.table(df, file= paste(outputDir,'/RNAseq/adonis_drugs_univariate_LiverRNA.txt',sep=""), sep="\t", quote=F,row.names = F)



