library(scales)
library(dplyr)

prs <- read.table("ad_from_allscore_best.finaltop14.txt", header=T)


rescale.many <- function(dat, column.nos, mix, max) { 
  nms <- names(dat) 
  for(col in column.nos) { 
    name <- paste(nms[col],".rescaled", sep = "") 
    dat[name] <- rescale(dat[,col], to=c(mix, max)) 
  } 
  cat(paste("Rescaled ", length(column.nos),      " variable(s)n")) 
  dat 
} 
id_diag <- subset(prs, select=c(1:2))
T <- subset(prs, select=-c(1:2))
l=length(names(prs))-2
rescaled <- rescale.many(T, c(1:l), -1,1)
rescaled <- subset(rescaled, select=-c(1:l))   
prs_rescaled=cbind(id_diag,rescaled)
#write.csv(prs_rescaled, "ADc12_microbiome_top14_rescaled.csv", row.names=F, quote=F)


#merge with info
dfinfo=read.csv("ADc12_IID_Diag_SEX_Age_APOE.csv")
dfmerge=left_join(dfinfo,prs_rescaled,by=c('IID','Diag'))
write.csv(dfmerge, "ADc12_microbiome_top14_withinfo_rescaled.csv", row.names=F, quote=F)