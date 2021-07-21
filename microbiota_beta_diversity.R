################################################################
# To-do:  beta diversity analysis of microbiota
# Author: Peishun Li, peishun@chalmers.se
################################################################
library(vegan)
library(dplyr)  #%>%
library(ggplot2)
library(extrafont)

# Set ggplot2 drawing parameter, such as axis line and text size, lengend and title size, and so on.
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
BARIA_wd <- "~/Documents/BARIA/baria"
outputDir <- "~/Documents/BARIA/baria/output"

################################
#  calculate the beta diversity
################################
#principal coordinates analysis PCoA - taxa/taxon/sample
comm <- speciePer
brayDistance<-vegdist(comm, method="bray", binary=FALSE) 
quantile(brayDistance)
brayDistance1<-as.matrix(brayDistance)
solPCoA <-cmdscale(brayDistance1,k=3,eig=TRUE)
eig <- solPCoA$eig
## ordination/ordinate of samples
ordi_site  <- data.frame(x = solPCoA$points[,1], y = solPCoA$points[,2], group=microbDesign106$group, sampleId=microbDesign106$sampleId)

# create a scatterplot with multi centroids (one for each class) that includes error bars
# The centroids should be positioned at the mean values for x and y for each class
# aggregate(ordi_site[,2:3],list(group=ordi_site$group),mean)
centroids <- aggregate(cbind(x,y)~group,ordi_site,mean)
# add centroids into ordi_site
ordi_site <- merge(ordi_site,aggregate(cbind(mean_x=x,mean_y=y)~group,ordi_site,mean),by="group")

## Data frame df_ell contains values to show ellipses. It is calculated with 
## function veganCovEllipse which is hidden in vegan package.
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) { 
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov))) }
#Generate ellipse points
df_ell <- data.frame()
for(g in levels(ordi_site$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(ordi_site[ordi_site$group==g,],
            veganCovEllipse(ord[[g]]$cov*sum(ordi_site$group==g),ord[[g]]$center))),group=g)) }
colnames(df_ell)<-c('x', 'y', 'group')

#######
## scatterplot has sample points, centroids and ellipses
p <- ggplot(data = ordi_site, aes(x, y,color =group)) + main_theme+ geom_point(alpha=.7, stroke = 0, size=0.9)+
#  geom_segment(aes(x=mean_x, y=mean_y, xend=x, yend=y), colour = "grey",size=0.2)+
  geom_point(data=centroids, stroke = 0,size=1.3,aes(shape=group))+
  geom_path(data=df_ell, aes(x=x, y=y,colour=group), size=0.3, linetype=2)+
labs(x=paste("PCoA 1 (", format(100 * eig[1] / sum(eig), digits=4), "%)", sep=""),
     y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""))+
  scale_colour_manual(values=c("#2c7bb6", "#fdae61", "#d7191c"))+
  scale_shape_manual(values=c( 17,17,17))+theme(legend.position="top")
ggsave(paste0(outputDir,"/microbiota/pcoa_species_bray_center.pdf"), p, width = 4, height = 3.5, units="cm", useDingbats=FALSE)


## PERMANOVA significance test for group-level differences
## Permutational multivariate analysis of variance (PERMANOVA) - Adonis (vegan)
## for community-level multivariate comparisons
set.seed(123)
## bray: sig, jaccard: sig, manhattan: no sig, euclidean: no sig,  chao: no sig
adonis_group <- adonis(comm ~ group, data = microbDesign106, permutations=999, method = "bray")
## or distance matrices
## adonis_group <- adonis(brayDistance~group, data=microbDesign106) 
print(as.data.frame(adonis_group$aov.tab)["group", "Pr(>F)"])
adonis_group$aov.tab$`Pr(>F)`[1]     ## 0.043 indicates significance



