################################################################
# Name: AUC.R
# To-do: calculate incremental AUC and total AUC of MMT
# Date: 27-8-2018
# Author: Peishun Li, peishun@chalmers.se
################################################################
AUC <- function(a,b){
  library(pracma)
  totalAUC <- trapz(a, b)
  iAUC <- vector()
  for (i in 1:(length(a)-1)){
    if((b[i+1]-b[1] >= 0) && (b[i]-b[1] >= 0))
    {iAUC[i] <-((b[i]-b[1]+b[i+1]-b[1])/2)*(a[i+1]-a[i])} 
    else if((b[i+1]-b[1] < 0) && (b[i]-b[1] >= 0))
    {iAUC[i] <-(b[i]-b[1])*((b[i]-b[1])/(b[i]-b[i+1])*(a[i+1]-a[i])/2)} 
    else if((b[i+1]-b[1] >= 0) && (b[i]-b[1] < 0))
    {iAUC[i] <-(b[i+1]-b[1])*((b[i+1]-b[1])/(b[i+1]-b[i])*(a[i+1]-a[i])/2)} 
    else if((b[i]-b[1] < 0) && (b[i+1]-b[1] < 0))
    {iAUC[i] <- 0}
  }
incrementalAUC<-sum(iAUC)
data.frame(totalAUC, incrementalAUC)
  
}