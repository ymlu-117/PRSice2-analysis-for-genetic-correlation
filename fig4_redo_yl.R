
#boxplot
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(EnvStats)
library(gridExtra)
library(grid)

setwd("E:/Rscript")

df=read.csv("ADc12_microbiome_top20_withinfo_Z_scored.csv")
df$rs429358_C=as.factor(df$rs429358_C)
#df$rs7412_T=as.factor(df$rs7412_T)

df$Diag_factor[df$Diag==1]<-"Case"
df$Diag_factor[df$Diag==0]<-"Control"
df$Diag_factor<-as.factor(df$Diag_factor)

dfid=df[,c(1:6,ncol(df))]
dfbest=df[,grepl( "best" , names(df) )]
dfbest=cbind(dfid,dfbest)

#rs429358_C
apoepos=read.table("apoe_pos_Zscore.txt",header=TRUE)
apoepos=apoepos$micro
dfapoepos=dfbest[,(names(dfbest) %in% apoepos)]
dfapoepos=cbind(dfid,dfapoepos)
colnames(dfapoepos)[8]<- "Collinsella_best"

head(dfapoepos)

yval=colnames(dfapoepos)[8]
yval_strip<-str_replace_all(yval,"_best","")
ymin<- min(dfapoepos$Collinsella_best)
ymax<- max(dfapoepos$Collinsella_best)
p=ggplot(dfapoepos, aes_string(x="rs429358_C", y=yval,fill="rs429358_C")) +
    geom_boxplot(outlier.color = NA)+
    xlab("APOE rs429358 Genotype")+
    ylim(ymin, ymax)+
    scale_x_discrete(labels=c("TT","TC","CC"))+
    ylab(paste("PRSs for",yval_strip,sep=" "))+
    stat_n_text(size=8,y.pos=-Inf,vjust=-0.5, fontface = "bold")+#n number
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())+
    theme(legend.position="none",
          axis.text.x = element_text(size = 20,face="bold",color="black"),#x axis label
          axis.title = element_text(size = 20,face="bold"),#axis title
          axis.text.y = element_text(size = 20,face="bold",color="black"))#y axis label
p=p + stat_compare_means(method = "anova",size=8)
outfilename=paste(yval_strip,"rs429358_C_boxplot.png",sep="_")
ggsave(outfilename,p,height=180,width=180,units="mm",dpi=300)

  