library(tidyverse)
df=read.csv("ADc12_IID_Diag_SEX_Age_APOE_Microbiometop14PRS77.csv")
dfid=df[,c(1,2)]
dfbest=df[,grepl( "best" , names(df) )]
dfbest=cbind(dfid,dfbest)

#can do this in the future:
#dfbest= select(df,IID,Diag,contains("best"))


dfout <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(dfout) <- c('coefficient', 'se', 'z-value','p','OR', '2.5%CI', '97.5.5%CI')

for (i in 3:ncol(dfbest)){
  modelglm=as.formula(paste0("Diag~",colnames(dfbest)[i]))
  glmfit<- glm(modelglm, data = dfbest, family = binomial)
  coefdf=as.data.frame(summary(glmfit)$coefficients)[2,]
  dfor=as.data.frame(exp(cbind(OR = coef(glmfit), confint(glmfit))))[2,]
  dfor_coef=cbind(coefdf,dfor)
  print(modelglm)
  dfout=rbind(dfout,dfor_coef)
}
colnames(dfout) <- c('coefficient', 'se', 'z-value','p','OR', '2.5%CI', '97.5.5%CI')
dfout=rownames_to_column(dfout,var="Traits")
dfout
write.table(dfout,file="microbiome_2ndtop14_regression06072022.csv",sep=",",quote=FALSE,row.names = FALSE)
