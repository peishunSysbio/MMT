##############################
# To-do: train predictive model using ridge regression
# Author: Peishun Li, peishun@chalmers.se
##############################
BARIA_wd <- "~/Documents/BARIA/baria"
outputDir <- "~/Documents/BARIA/baria/output"

library(glmnet)
library(ggpubr)

## using all data: X_data, Y_data
set.seed(n=1)
cv_fit <- cv.glmnet(X_data, Y_data, alpha = 0, nfolds = 10)
# cv_fit <- cv.glmnet(datX, datY, nlambda = 20, alpha = 0, nfolds = 10, intercept = int, standardize.response = sr)
best_lamda <- cv_fit$lambda.min

# build the model
# cv_fit$glmnet.fit
best_ridge <- glmnet(X_data, Y_data, alpha = 0, lambda = best_lamda)
# summary(best_ridge)

######## To evaluate performance of the ridge model with the optimal lambda, 
## 5-fold cross-validation (samples were randomly divided into five equal parts. 
kfoldCV <- function(datX, datY, k, a = 0, int = TRUE, sr = FALSE, n) {
  set.seed(n)
  kgrp <- split(1:length(datY), sample(1:k, length(datY), replace = T))
  ## table(sample(1:k, length(datY), replace = T))
  rmse <- array(dim = k)
  scor <- rmse
  pcor <- rmse
  for (i in 1:k) {
    test_i <- unlist(kgrp[[i]])
    fit_i <- glmnet(datX[-test_i, ], datY[-test_i], lambda = best_lamda, alpha = a, intercept = int, standardize.response = sr)
    predict_i <- predict(fit_i, datX[test_i, ], s = best_lamda, alpha = a, intercept = int, standardize.response = sr)
    if (sr == TRUE) {
      rmse[i] <- sqrt(sum((predict_i - scale(datY[test_i]))^2) / length(predict_i))
      scor[i] <- cor(predict_i, scale(datY[test_i]), method = "spearman")
    } else {
      rmse[i] <- sqrt(sum((predict_i - datY[test_i])^2) / length(predict_i))
      temp <- cor.test(predict_i, datY[test_i], method = "spearman")
      scor[i] <- as.numeric(temp$estimate)
      pcor[i] <- as.numeric(temp$p.value)
    }
  }
  return(cbind(rmse, scor, pcor))
}
#
resKfoldCV <- kfoldCV(X_data, Y_data, k = 5, a = 0, n=n)

## met0h
met0h_ridge_regression_model <- resKfoldCV 

## met2h
met2h_ridge_regression_model <- resKfoldCV 

##
df <- data.frame(RMSE=resKfoldCV[,1], type=rep('metR',5))
RMSE_plot<- ggplot(df, aes(x=type, y=RMSE, colour =type)) + geom_boxplot(width=0.4, outlier.shape = NA) +
  labs(x='',y="RMSE", title="")+ main_theme +theme(legend.position = 'none' ) +
  geom_jitter(width = 0.2,size=0.5)
#
ggsave(RMSE_plot, filename = paste0(data_dir, "/rmse_5foldCV_ridge_predict.pdf"), width = 4, height = 4, units = "cm")


#######  predict using the all data.
model_predict <- predict(best_ridge, s = best_lamda, newx = X_data)
#
prediction_pairs <- data.frame(original = as.vector(Y_data),
                               predict = as.vector(model_predict))

outFigure <- paste0(data_dir, "/cv_fit_train_met0h_ridge_predict.pdf")
pdf(outFigure, width=4, height=4,useDingbats=FALSE)
plot(cv_fit, main = paste0("set.seed(", n, "), best_lamda=",best_lamda))
dev.off()


# color="#d7301f", pearson, spearman
prediction_plot <- ggplot(prediction_pairs, aes(x=original, y=predict)) +
  geom_point(size=0.5, color="black") +
  geom_smooth(method=lm, se = T) + stat_cor(method = "spearman", label.y = 1100, size=2)+  #label.x = 300,
  labs(x='Actual glucose tAUC',y='Predictive glucose tAUC', title="")+ main_theme
#
ggsave(prediction_plot, filename = paste0(data_dir, "/accuracy_train_met0h_ridge_predict_lm.pdf"), width = 4, height = 4, units = "cm", useDingbats=FALSE)

##### feature importance
# predict(best_ridge,type="coef")
best_ridge_coef <- coef(best_ridge)
best_ridge_coef <- as.matrix(coef(best_ridge))
best_ridge_coef <- data.frame(coef = best_ridge_coef[,1], name = rownames(best_ridge_coef))[-1,]
summary(best_ridge_coef$coef)
best_ridge_coef$coef[abs(best_ridge_coef$coef)>=1]
#
best_ridge_coef <- merge(best_ridge_coef, metabAnno[,2:4],  by='row.names', all.x = T)
best_ridge_coef_sort <- best_ridge_coef[rev(order(abs(best_ridge_coef$coef))),]
#
write.table(best_ridge_coef_sort, file= paste(outputDir,'/ML_regression/met0hCoef.txt',sep=""), sep="\t", quote=F,row.names = F)

## Show coefficients for the top 30 mets at 0h
top30met0hCoef <- best_ridge_coef[rev(order(abs(best_ridge_coef$coef))),][1:32, ]
top30met0hCoef <- top30met0hCoef[!top30met0hCoef$SUB.PATHWAY=='',]
top30met0hCoef <- top30met0hCoef[order(top30met0hCoef$coef),]

## fasting
df <- data.frame(mets=top30met0hCoef$BIOCHEMICAL, value=top30met0hCoef$coef,Class =top30met0hCoef$SUPER.PATHWAY)
mycols <- c("#b3e2cd", "#fdcdac", "#cbd5e8", "#e6f5c9", "#fff2ae", "#f1e2cc","#F39B7F99")
df$Class <- factor(df$Class , levels = c('Lipid', 'Carbohydrate', 'Amino Acid','Nucleotide', 'Xenobiotics',  'Peptide','Energy'))

##
df$mets<-factor(df$mets,levels = df$mets)
p <- ggplot(data = df)+geom_bar(aes(x=mets,y=value, fill=Class), stat="identity", colour='white', size=0.1,alpha = 1)+
  coord_flip()+main_theme+labs(y='coefficients')+ scale_fill_manual(values = mycols) +  
theme(legend.position=c(0.9,0.2), axis.text.x =element_text(size=5), axis.text.y=element_text(size=5), 
      legend.text= element_text(size=5), text=element_text(family="Arial", size=6), legend.key.size = unit(0.2, "cm"))
##
ggsave(paste(data_dir, "/best_ridge_met0h_top30Coef.pdf", sep=""), p, width = 7, height =7, units = "cm")


