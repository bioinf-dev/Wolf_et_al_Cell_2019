
library(plyr)
require(reshape2)
library(gridExtra)
library(cowplot)
library(survival)
library('survminer')
library(scales)
library('dplyr')
rm(list=ls(all=TRUE))

#Snyder data
dat<-read.table("yardena/Riaz_analysis/final_data/SDI_clinical_data_Snyder.txt",header=T,sep="\t")
timestrata.surv1 <- survfit(Surv(time,censor)~ SDI_group, dat)
km1<-ggsurvplot(timestrata.surv1,risk.table = TRUE, risk.table.height=0.25,risk.table.title="",font.x=9,font.y=9,font.tickslab=10,legend = "none",xlab="Time (months)",ylab="Progression Free Survival (%)",tables.theme = theme_cleantable(),palette = c("#D55E00","#009E73"))
summary(coxph(formula=Surv(time, censor)~ SDI, dat))

#Riaz data
dat<-read.table("yardena/Riaz_analysis/final_data/SDI_clinical_data_Riaz.txt",header=T,sep="\t")
timestrata.surv1 <- survfit(Surv(time,censor)~ SDI_group, dat)
km2<-ggsurvplot(timestrata.surv1,risk.table = TRUE, risk.table.height=0.25,risk.table.title="",font.x=9,font.y=9,font.tickslab=10,legend = "none",xlab="Time (months)",ylab="Progression Free Survival (%)",tables.theme = theme_cleantable(),palette = c("#D55E00","#009E73"))
summary(coxph(formula=Surv(time, censor)~ SDI, dat))

#Hugo data
dat<-read.table("yardena/Riaz_analysis/final_data/SDI_clinical_data_Hugo.txt",header=T,sep="\t")
timestrata.surv1 <- survfit(Surv(time,censor)~ SDI_group, dat)
km3<-ggsurvplot(timestrata.surv1,risk.table = TRUE, risk.table.height=0.25,risk.table.title="",font.x=9,font.y=9,font.tickslab=10,legend = "none",xlab="Time (months)",ylab="Progression Free Survival (%)",tables.theme = theme_cleantable(),palette = c("#D55E00","#009E73"))
summary(coxph(formula=Surv(time, censor)~ SDI, dat))

#Van Allen data
dat<-read.table("yardena/Riaz_analysis/final_data/SDI_clinical_data_VanAllen.txt",header=T,sep="\t")
timestrata.surv1 <- survfit(Surv(time,censor)~ SDI_group, dat)
km4<-ggsurvplot(timestrata.surv1,risk.table = TRUE, risk.table.height=0.25,risk.table.title="",font.x=9,font.y=9,font.tickslab=10,legend = "none",xlab="Time (months)",ylab="Progression Free Survival (%)",tables.theme = theme_cleantable(),palette = c("#D55E00","#009E73"))
summary(coxph(formula=Surv(time, censor)~ SDI, dat))


splots <- list()
splots[[1]]<-km1 
splots[[2]]<-km2
splots[[3]]<-km3
splots[[4]]<-km4
arrange_ggsurvplots(splots, print = TRUE,ncol = 2,nrow=2)

# Forest plot
ORs<-read.table("yardena/Riaz_analysis/final_data/for_forest.txt",header=TRUE,sep='\t')
ORs$Study <- factor(ORs$Study, levels=(ORs$Study))

fp <- ggplot(data=ORs,aes(x=Study, y=OR, ymin=lower, ymax=upper))+geom_pointrange(colour="darkblue",size=1.0) + geom_hline(yintercept=1, lty=2)+scale_y_log10(breaks=c(0.1,1.0,10.0),limits=c(0.1,100)) +xlab("Risk events") + ylab("Hazard Ratio (95% CI) (log scale)")+ theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, hjust = 1))

fp+theme_classic()+coord_flip()

library("CombinePValue")
selfcontained.test(pvalue=c(0.0064,0.079,0.096,0.96),weight=NA,p_permu=NA)


