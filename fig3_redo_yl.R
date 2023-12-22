
library(ggpubr)
library(tidyverse)
library(ggplot2)
library(EnvStats)
library(gridExtra)
library(grid)

#change here: to your working directory. You'll need to input dir like shown below.
setwd("H:/dbGaP_data/bioinfo_guides/microbiome_proj_script")
#use command below to make sure your working dir is right.
                          
df=read.csv("ADc12_microbiome_top20_withinfo_Z_scored.csv")#change input file name if necessary

df$Diag_factor[df$Diag==1]<-"Case"
df$Diag_factor[df$Diag==0]<-"Control"
df$Diag_factor<-as.factor(df$Diag_factor)

dfid=df[,c(1:6,ncol(df))]
dfbest=df[,grepl( "best" , names(df) )]
dfbest=cbind(dfid,dfbest)

pos=read.table("casecontrol_pos_Zscore.txt",header=TRUE)
neg=read.table("casecontrol_neg_Zscore.txt",header=TRUE)

pos=pos$micro
neg=neg$micro

dfpos=dfbest[,pos]
dfneg=dfbest[,neg]
dfpos=cbind(dfid,dfpos)
dfneg=cbind(dfid,dfneg)

#neg group


test=dfneg
wilcox=data.frame(micro=names(test)[8:(ncol(test))],p=NA)
for (i in 8:(ncol(test))){
  nr=i-7
  yval=colnames(test)[i]
  yval_strip<-str_replace_all(yval,"_best","")
  modelfit=as.formula(paste0(colnames(test)[i]," ~ Diag_factor"))
  res=wilcox.test<-wilcox.test(modelfit, data = test,exact=FALSE)
  pval=res$p.val
  wilcox$p[nr]=pval
}
wilcox$p<-signif(wilcox$p,digits=3)

wilcox$p

ADcase_plot <- function(i,test,ymin,ymax) {
  yval=colnames(test)[i]
  yval_strip<-str_replace_all(yval,"_best","")
  if (grepl("Eubacterium",yval_strip)) {
    yval_strip=str_replace_all(yval_strip,"Eubacterium","E.")
    yval_strip=str_replace_all(yval_strip,"group"," group")
  }
  pval= wilcox$p[wilcox$micro==yval]
  pval=formatC(pval,format="e",digits=2)
  pval_str=paste0("p=",pval)
  p=ggplot(test, aes_string(x="Diag_factor", y=yval,fill="Diag_factor")) +
    geom_boxplot(outlier.color = NA)+ 
    annotate("text",  label= pval_str, x=-Inf, y = Inf, vjust=1.2, hjust=-0.1, size=5,fontface ="bold") +#pval text
    ylim(-3.0,3.0)+
    labs(title=yval_strip)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()
          ,plot.title = element_text(hjust = 0.5,size=15,face="bold"))+#title
    theme(legend.position="none",
          axis.text.y = element_text(size = 15,face="bold"), #y axis
          axis.text.x = element_text(size = 17,face="bold",color="black"), #x axis
          axis.title = element_blank())
  return(p)
}

#Adlercreutzia_best
#run  i from 8 to 13
i=8 # change here Adlercreutzia  
n=i-7
#ymin<- min(dfneg$Adlercreutzia_best)
#ymax<- max(dfneg$Adlercreutzia_best)
ymin<- -3.0
ymax<- 3.0
p=ADcase_plot(i,test,ymin,ymax)
assign(paste0("neg_p",n),p)
neg_p1

#Eisenbergiella_best
i=9 
n=i-7 #change here Eisenbergiella_best
#ymin<- min(dfneg$Eisenbergiella_best)
#ymax<- max(dfneg$Eisenbergiella_best)
ymin<- -3.0
ymax<- 3.0
p=ADcase_plot(i,test,ymin,ymax)
assign(paste0("neg_p",n),p)
neg_p2

#Eubacteriumfissicatenagroup
i=10 
n=i-7 #change here Eubacteriumfissicatenagroup  
#ymin<- min(dfneg$Eubacteriumfissicatenagroup_best)
#ymax<- max(dfneg$Eubacteriumfissicatenagroup_best)
ymin<- -3.0
ymax<- 3.0
p=ADcase_plot(i,test,ymin,ymax)
assign(paste0("neg_p",n),p)
neg_p3
i=11 
n=i-7 #change here Eubacteriumnodatumgroup_best
ymin<- min(dfneg$Eubacteriumnodatumgroup_best)
ymax<- max(dfneg$Eubacteriumnodatumgroup_best)
p=ADcase_plot(i,test,ymin,ymax)
assign(paste0("neg_p",n),p)
neg_p4
i=12 
n=i-7 #change here Gordonibacter  
ymin<- min(dfneg$Gordonibacter_best)
ymax<- max(dfneg$Gordonibacter_best)
p=ADcase_plot(i,test,ymin,ymax)
assign(paste0("neg_p",n),p)
neg_p5

i=13 
n=i-7 #change here Prevotella9  
#ymin<- min(dfneg$Prevotella9_best)
#ymax<- max(dfneg$Prevotella9_best)
ymin<- -3.0
ymax<- 3.0
p=ADcase_plot(i,test,ymin,ymax)
assign(paste0("neg_p",n),p)
neg_p6


#after the run you should have neg_p1 to neg_p6
a=grid.arrange(neg_p1,neg_p2,neg_p3,neg_p4,neg_p5,neg_p6,ncol=3)
a
ggsave("ADcase_neg_Zscore_New.png",a,height=150,width=180,units="mm",dpi=300)



#-------------------------

test=dfpos

wilcox=data.frame(micro=names(test)[8:(ncol(test))],p=NA)
for (i in 8:(ncol(test))){
  nr=i-7
  yval=colnames(test)[i]
  yval_strip<-str_replace_all(yval,"_best","")
  modelfit=as.formula(paste0(colnames(test)[i]," ~ Diag_factor"))
  res=wilcox.test<-wilcox.test(modelfit, data = test,exact=FALSE)
  pval=res$p.val
  wilcox$p[nr]=pval
}
wilcox$p<-signif(wilcox$p,digits=3)
ADcase_plot <- function(i,test,ymin,ymax) {
  yval=colnames(test)[i]
  yval_strip<-str_replace_all(yval,"_best","")
  pval= wilcox$p[wilcox$micro==yval]
  pval=formatC(pval,format="e",digits=2)
  pval_str=paste0("p=",pval)
  p=ggplot(test, aes_string(x="Diag_factor", y=yval,fill="Diag_factor")) +
    geom_boxplot(outlier.color = NA)+ 
    annotate("text",  label= pval_str, x=-Inf, y = Inf, vjust=1.2, hjust=-0.1, size=5,fontface ="bold") +#pval text
    ylim(-3.0,3.0)+
    labs(title=yval_strip)+
    theme_bw()+
    theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()
          ,plot.title = element_text(hjust = 0.5,size=15,face="bold"))+#title
    theme(legend.position="none",
          axis.text.y = element_text(size = 15,face="bold"), #y axis
          axis.text.x = element_text(size = 17,face="bold",color="black"), #x axis
          axis.title = element_blank())
  return(p)
}



#risk
#Bacteroides_best
#Collinsella_best
#Lachnospira_best
#Veillonella_best


#pos group for each
#run i from 8 to 13, now change: 8 to 11 

i=8#change here for Bacteroides
n=i-7
#ymin<- min(dfpos$Bacteroides_best)
#ymax<- max(dfpos$Bacteroides_best)
ymin<- -3.0
ymax<- 3.0
p=ADcase_plot(i,test,ymin,ymax)
assign(paste0("pos_p",n),p)
pos_p1


i=9#change here for Collinsella
n=i-7
#ymin<- min(dfpos$Collinsella_best)
#ymax<- max(dfpos$Collinsella_best)
ymin<- -3.0
ymax<- 3.0
p=ADcase_plot(i,test,ymin,ymax)
assign(paste0("pos_p",n),p)
pos_p2

i=10#change here for Lachnospira_best
n=i-7

#ymin<- min(dfpos$Lachnospira_best)
#ymax<- max(dfpos$Lachnospira_best)
ymin<- -3.0
ymax<- 3.0
p=ADcase_plot(i,test,ymin,ymax)
assign(paste0("pos_p",n),p)
pos_p3


i=11#change here for Veillonella
n=i-7
#ymin<- min(dfpos$Veillonella_best)
#ymax<- max(dfpos$Veillonella_best)
ymin<- -3.0
ymax<- 3.0
p=ADcase_plot(i,test,ymin,ymax)
assign(paste0("pos_p",n),p)
pos_p4


#after the run you should have pos_p1 to pos_p4
a=grid.arrange(pos_p1,pos_p2,pos_p3, pos_p4, ncol=2)
a
ggsave("ADcase_pos_Zscore_New.png",a,height=180,width=180,units="mm",dpi=300)



