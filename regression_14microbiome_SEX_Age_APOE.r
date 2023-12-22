library(tidyverse)
df=read.csv("ADc12_IID_Diag_SEX_Age_APOE_Microbiometop14PRS77.csv")
dfid=df[,c(1:6)]
dfbest=df[,grepl( "best" , names(df) )]
dfbest=cbind(dfid,dfbest)

#can do this in the future:
#dfbest= select(df,IID,Diag,contains("best"))

dfout <- data.frame(matrix(ncol = 8, nrow = 0))
colnames(dfout) <- c('Trait','coefficient', 'se', 'z-value','p','OR', '2.5%CI', '97.5.5%CI')

for (i in 7:ncol(dfbest)){
 modelglm=as.formula(paste0("Diag~", "SEX.rescaled + Age.rescaled + rs429358_C.rescaled + rs7412_T.rescaled + ",  colnames(dfbest)[i]))
	
  micro=str_replace(colnames(dfbest)[i],"_best.rescaled","")
  glmfit<- glm(modelglm, data = dfbest, family = binomial)

  coefdf=as.data.frame(summary(glmfit)$coefficients)[-1,]
  dfor=as.data.frame(exp(cbind(OR = coef(glmfit), confint(glmfit))))[-1,]
  dfor_coef=cbind(coefdf,dfor)
  dfor_coef=rownames_to_column(dfor_coef,var='Trait')
  for (n in 1:nrow(dfor_coef)-1){
	dfor_coef[n,1]=paste(micro,dfor_coef[n,1],sep="_")
  }
  print(modelglm)
  dfout=rbind(dfout,dfor_coef)
}
colnames(dfout) <- c('coefficient', 'se', 'z-value','p','OR', '2.5%CI', '97.5.5%CI')

write.table(dfout,file="microbiome_top14_SEX_Age_APOEregression06072022.csv",sep=",",quote=FALSE,row.names = FALSE)

#model: prs vs 2 apoe vaviants
