################################################################
# To-do: analyse mixed meal test (MMT) data
# Author: Peishun Li, peishun@chalmers.se
################################################################
library(ggplot2)
library(reshape)
library(ggpubr)
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
###################################################################
BARIA_wd <- "~/Documents/BARIA/baria"
outputDir <- "~/Documents/BARIA/baria/output"

###################################################################
## All MMT vars of 106 patients with multi-omics
mmtAllVars106 <- mmtAllVars[rownames(metadata106),]

# import MMT-related vars including Apo A1, Apo B, Chol,Tri,HDLc	Lp(a),NEFA from 157 patient ids
mmtLipiRaw <- read.table(paste0(BARIA_wd,"/data/clinicalDB/triglycerides.txt"), header=T, sep="\t",stringsAsFactors=F, dec='.')
#
mmtLipidata <- mmtLipiRaw[,c(1,3:5,7:9)]
mmtLipidata <- mmtLipidata[!is.na(mmtLipidata$Pt.ID),]
rownames(mmtLipidata) <- mmtLipidata$Pt.ID
mmtLipidata106 <- mmtLipidata[rownames(mmtAllVars106),]
mmtLipidata106 <- mmtLipidata106[,-1]
mmtLipidata106$Chol <-as.numeric(mmtLipidata106$Chol)
mmtLipidata106$Lp.a. <-as.numeric(mmtLipidata106$Lp.a.)
# Triglycerides (mmol/L)
trigMMTdata <- mmtLipiRaw[c(1,2,6)]
trigMMTdata$id <- rep(mmtLipidata$Pt.ID,each=4)
trigMMTdata <- trigMMTdata[,-1]
trigMMTdata1 <- reshape(trigMMTdata, idvar = "id", timevar = "t", direction = "wide")
rownames(trigMMTdata1) <- trigMMTdata1$id
trigMMTdata106 <- trigMMTdata1[rownames(mmtAllVars106),]
trigMMTdata106 <-trigMMTdata106[,-1]

#
identical(rownames(medicationInf106), rownames(mmtLipidata106))
df <- mmtLipidata106
df <- cbind(trigMMTdata106,allAUC_trig)
df$group <-medicationInf106$groupFinal
df_long <- melt(df, id='group')
##
p <- ggplot(df_long, aes(x=group, y=value, color=group)) +
  geom_boxplot(alpha=0.4,outlier.size = 0) +facet_wrap(~variable,scales ="free", ncol = 3)+
  stat_summary(fun.y=mean, geom="point", shape=3, size=1) +
  main_theme+ scale_colour_manual(values=c("#2c7bb6", "#fdae61", "#d7191c"))+
  stat_compare_means(method = "kruskal.test",label = "p.format", label.x = 2,size = 1.5)+
  stat_compare_means(comparisons = my_comparisons,method ='wilcox.test',label = "p.signif",size = 1.5)+
  labs(x='',y='Level', title="")+theme(legend.position="non")
# show_point_shapes(), outlier.shape = NA
ggsave(paste(outputDir,'/clinic/MMTT_trigAUC_Boxplot.pdf',sep=""), p, width = 4, height = 4)

################### impute missing data for calculate AUC
library(DMwR)
glucMMTImputation<-knnImputation(mmtAllVars106[,1:7])
insulMMTImputation<-knnImputation(mmtAllVars106[,8:14])
trigMMTImputation<-knnImputation(trigMMTdata106)

## Area under the curve (AUC)
source("AUC.R")
########################### AUC of 106 patients with MMT
allAUC <- data.frame()
a <- c(0, 10, 20, 30, 60, 90, 120)
for(j in 1:nrow(mmtAllVars106)) {
  #b<-as.numeric(glucMMTImputation[j,])    ## glucose
  b<-as.numeric(insulMMTImputation[j,])   ## insulin
  eachVar <- AUC(a,b)
  rownames(eachVar) <- j
  allAUC<- rbind(allAUC,eachVar)
  #
}
allAUC_gluc <- allAUC
colnames(allAUC_gluc)<-c("totalAUC_gluc","incrementalAUC_gluc")
allAUC_insulin <- allAUC
colnames(allAUC_insulin)<-c("totalAUC_insulin","incrementalAUC_insulin")
# for Triglycerides
allAUC <- data.frame()
a <- c(0, 30, 60, 120)
for(j in 1:nrow(mmtAllVars106)) {
  b<-as.numeric(trigMMTImputation[j,])
  eachVar <- AUC(a,b)
  rownames(eachVar) <- j
  allAUC<- rbind(allAUC,eachVar)
}
allAUC_trig <- allAUC
colnames(allAUC_trig)<-c("totalAUC_trig","incrementalAUC_trig")
#
mmtAllVars106 <- cbind(mmtAllVars106, cbind(allAUC_gluc, allAUC_insulin))

# Estimates of beta-cell function
# Insulinogenic index (pmol/mmol, an estimate of early insulin secretion) was calculated by dividing the increment 
# in insulin during the first 30 postprandial minutes by the increment in glucose over the same period Insulinogenic_index
allAUC30 <- data.frame()
a <- c(0, 10, 20, 30)
for(j in 1:nrow(mmtAllVars106)) {
  #b<-as.numeric(glucMMTImputation[j,c(1:4)])    
  b<-as.numeric(insulMMTImputation[j,c(1:4)])  
  eachVar <- AUC(a,b)
  rownames(eachVar) <- j
  allAUC30<- rbind(allAUC30,eachVar)
}
allAUC30_gluc <- allAUC30
colnames(allAUC30_gluc)<-c("totalAUC_gluc30","incrementalAUC_gluc30")

allAUC30_insulin <- allAUC30
colnames(allAUC30_insulin)<-c("totalAUC_insulin30","incrementalAUC_insulin30")
mmtAllVars106 <- cbind(mmtAllVars106, cbind(allAUC30_gluc,allAUC30_insulin))
##
mmtAllVars106[mmtAllVars106$incrementalAUC_gluc30==0,'incrementalAUC_gluc30']=1
mmtAllVars106$insulinogenic_index <- mmtAllVars106$incrementalAUC_insulin30/mmtAllVars106$incrementalAUC_gluc30

# AUCinsulin/AUCglucose ratio (pmol/mmol)
mmtAllVars106$AUCinsulin_AUCglc <- mmtAllVars106$totalAUC_insulin/mmtAllVars106$totalAUC_gluc
# iAUCinsulin/iAUCglucose (pmol/mmol)
mmtAllVars106$iAUCinsulin_iAUCglc <- mmtAllVars106$incrementalAUC_insulin/mmtAllVars106$incrementalAUC_gluc

# Estimates of peak time
a <- c(0, 10, 20, 30, 60, 90, 120)
mmtAllVars106$peak_time <- apply(glucMMTImputation,1,function(x) a[which.max(x)])
##
write.table(mmtAllVars106, file= paste0(BARIA_wd,"/data/clinicalDB/20190625/mmtAllVars106.txt"), sep="\t", quote=F)

############################################### Plot related to MMT of 106 patients 
## Plot glucose level in MMT of 106 people in multi groups
identical(row.names(mmtAllVars106),row.names(medicationInf106))
glucMealTest <- mmtAllVars106[,c(1:7)]
insulinMealTest <- mmtAllVars106[,c(8:14)]
MMTData<- glucMealTest
MMTData<- insulinMealTest
# Triglycerides (mmol/L)
summary(unlist(trigMMTdata106))
MMTData<- trigMMTdata106

##
MMTData$id <- rownames(MMTData)
MMTData$group <- medicationInf106$groupFinal
MMTDatae_long <- melt(MMTData, id=c("id","group"))
# for glucose and insulin (pmol/L)
MMTDatae_long$time <- rep(c(0,10,20,30,60,90,120), each=nrow(mmtAllVars106))
# for triglycerides mmol/l
MMTDatae_long$time <- rep(c(0,30,60,120), each=nrow(mmtAllVars106))
sum(is.na(MMTDatae_long$value))

# summarySE provides the standard deviation, standard error of the mean, and a (default 95%) confidence interval
source('/Users/peishun/Documents/BARIA/baria/src/summarySE.R')
MMTDataSE <- summarySE(MMTDatae_long, measurevar="value", groupvars=c("group","time"),na.rm=T)
# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(7) # move them .05 to the left and right
p <- ggplot(MMTDataSE, aes(x=time, y=value, colour=group)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=10, position=pd) +
  geom_line(position=pd) +geom_point(position=pd,size = .6)+ scale_colour_manual(values=c("#2c7bb6", "#fdae61", "#d7191c"))+
  main_theme+labs(x='Time (minutes)',y='Insulin level (pmol/l)', title="")+theme(legend.position="none")+
  scale_x_continuous(breaks=seq(0,120,30))
ggsave(paste0(outputDir,"/clinic/MMTGlucMeanse3g2.pdf"), p, width =5.5, height = 5, units = 'cm', useDingbats=FALSE)
ggsave(paste0(outputDir,"/clinic/MMTInsulinMeanse3g.pdf"), p, width =5.5, height = 5, units = 'cm', useDingbats=FALSE)  
ggsave(paste0(outputDir,"/clinic/MMTrigMeanse3g.pdf"), p, width =5.5, height = 5, units = 'cm', useDingbats=FALSE)  

# Plot for each patients : Triglycerides
p<- ggplot(data=MMTDatae_long, aes(x=time, y=value,group=id,colour=group)) + 
  geom_line(aes(colour=group),size=.2,alpha=0.9)+geom_point(size = 0.2,alpha=0.4)+
  scale_colour_manual(values=c("#2c7bb6", "#fdae61", "#d7191c")) +
  main_theme+labs(x='Time (min)',y='Triglycerides level (mmol/l)', title="")+
  theme(legend.position="none")+ scale_x_continuous(breaks=seq(0,120,30))
## geom_hline(yintercept = c(5.6,5.7, 6.5,7),color="black",linetype=5,size=0.4)
ggsave(paste0(outputDir,"/clinic/MMTGlucAll_1.pdf"), p, width =5, height = 5, units = 'cm', useDingbats=FALSE)
ggsave(paste0(outputDir,"/clinic/MMTInsulincAll_1.pdf"), p, width =5, height = 5, units = 'cm', useDingbats=FALSE)
ggsave(paste0(outputDir,"/clinic/MMTTrigAll_1.pdf"), p, width =5, height = 5, units = 'cm', useDingbats=FALSE)


## total insulin AUC vs glucose AUC
### Avoid overlapping labels in ggplot2 charts
install.packages('ggrepel')
library(ggrepel)
# geom_text_repel, geom_label_repel 
# spearman
#cor(df$totalAUC_gluc,df$totalAUC_insulin, method = "pearson")
p <- ggplot(df, aes(x=totalAUC_gluc, y=totalAUC_insulin/1000, color=group))+ geom_point(alpha=.7, size=0.7) +
  labs(x='tAUC of glucose (mmol/l*min)',y='tAUC of insulin (nmol/l*min)', title="")+
  #  geom_text_repel(aes(label=rownames(df)), size=1, segment.size = 0.2,segment.alpha =0.4) + 
  main_theme+ theme(legend.position="none")+geom_smooth(method = "lm", se = F)+
  stat_cor(method = "spearman", label.x = 1300,size=2)+
  scale_colour_manual(values=c("#2c7bb6", "#fdae61", "#d7191c"))
# 
ggsave(paste(outputDir,'/clinic/MMTT_TotalAUC_insulin_glucose-spearman.pdf',sep=""), p, width =5.5, height = 5, units = 'cm', useDingbats=FALSE)


## incremental insulin AUC vs glucose AUC --spearman
p = ggplot(df, aes(x=incrementalAUC_gluc, y=incrementalAUC_insulin/1000, color=group))+ geom_point(alpha=.7, size=0.7) +
  labs(x='iAUC of glucose (mmol/l*min)',y='iAUC of insulin (nmol/l*min)', title="")+ main_theme+ theme(legend.position="none")+
  scale_colour_manual(values=c("#2c7bb6", "#fdae61", "#d7191c"))+geom_smooth(method = "lm", se = FALSE)+
  stat_cor(method = "spearman",label.x = 300,size=2)
# title="iAUC insulin vs glucose")
ggsave(paste(outputDir,'/clinic/MMTT_iAUC_insulin_glucose-spearman.pdf',sep=""), p, width =5.5, height = 5, units = 'cm', useDingbats=FALSE)



#### anova of MMT data
df<-mmtAllVars106[,1:7] # glucose
df<-mmtAllVars106[,8:14] #insuline
df<-trigMMTdata106
#
identical(rownames(df), rownames(medicationInf106))
df$group <- factor(medicationInf106$groupFinal)
df$id <- factor(paste0('id',rownames(df)))
df_long <- melt(df, id=c("id","group"))

#
summary(aov(value ~  group * variable, data=df_long))
aov.out <-  summary(aov(value ~  group * variable + Error(id/variable), data=df_long))
# glucose: group p < 1.1e-15
# insulin time: p < 1.1e-15
# trigl: group p >0.05


 