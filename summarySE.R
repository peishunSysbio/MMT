################################################################
# To-do: Gives count, mean, standard deviation
# Author: Peishun Li, peishun@chalmers.se
################################################################
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    library(plyr)

    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
        if (na.rm) sum(!is.na(x))
        else       length(x)
    }

    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    # The Interquartile Range: IQR(x) = quantile(x, 3/4) - quantile(x, 1/4)
    datac <- ddply(data, groupvars, .drop=.drop,
      .fun = function(xx, col) {
        c(N    = length2(xx[[col]], na.rm=na.rm),
          mean = mean(xx[[col]], na.rm=na.rm),
          sd   = sd(xx[[col]], na.rm=na.rm),
          median   = median(xx[[col]], na.rm=na.rm),
          max   = max(xx[[col]], na.rm=na.rm),
          min   = min(xx[[col]], na.rm=na.rm),
          normTest= shapiro.test(xx[[col]])$p.value>0.05,
          downIQR = quantile(xx[[col]],probs =1/4, na.rm=na.rm),
          upIQR = quantile(xx[[col]],probs =3/4, na.rm=na.rm)
        )
      },
      measurevar
    )

    # Rename the "mean" column    
    datac <- rename(datac, c("mean" = measurevar))

    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval: 
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult

    return(datac)
}