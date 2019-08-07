library(GGally)
library(survival)
library(ggplot2)
library(ggpubr)
library(readr)
#Read in processed TCGA data related to figure 1
survdat <- read_delim("./Figure1Data.txt",delim = "\t", escape_double = FALSE, trim_ws = TRUE)
survdat = as.data.frame(survdat)
rownames(survdat) = survdat$Patient
survdat = survdat[,-1]


#TCGA survival analysis

surv = survfit(Surv(time,status) ~ MUTATION_LOAD, data = survdat)
diff = survdiff(Surv(time,status) ~ MUTATION_LOAD, data = survdat)

p2 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 2),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_manual(name = "Mutation load", values = c("blue", "red"), labels = c("blue" = "low","red" =  "high")) + scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p2)

pdf("./Fig1A.pdf")
print(p2)
dev.off()

cphMUTATION_LOAD = coxph(Surv(time,status) ~ age + stage + mutation_load, data = survdat[!is.na(survdat$mutation_load),])

surv = survfit(Surv(time,status) ~ CNV_LOAD, data = survdat)
diff = survdiff(Surv(time,status) ~ CNV_LOAD, data = survdat)

p2 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 2),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_manual(name = "CNV load", values = c("blue", "red"), labels = c("blue" = "low","red" =  "high")) + scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p2)

pdf("./Fig1B.pdf")
print(p2)
dev.off()

cphCNV_LOAD = coxph(Surv(time,status) ~ age + stage + cnv_load, data = survdat)


surv = survfit(Surv(time,status) ~ ITH, data = survdat)
diff = survdiff(Surv(time,status) ~ ITH, data = survdat)

p1 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 2),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_manual(name = "ITH(#clones)", values = c("blue", "red"), labels = c("blue" = "low","red" = "high")) + scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p1)
pdf("./Fig1C.pdf")
print(p1)
dev.off()

cphITH = coxph(Surv(time,status) ~ age + stage + purity + clones, data = survdat)


cols = c("#42a7f4", "#3863d8", "#31f3f9","#9e32fc","#ff30e6","#fc1919")
vals = c("2", "3","1","4","5","6")
surv = survfit(Surv(time,status) ~ clones, data = survdat)
diff = survdiff(Surv(time,status) ~ clones, data = survdat)
p1 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 3),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_manual(name = "#clones", values = cols, labels = vals) + scale_y_continuous(labels = scales::percent) + scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")

print(p1)
pdf("./Fig1D.pdf")
print(p1)
dev.off()


survdat$ITHvsMutationLoad = strata(survdat$ITH,survdat$MUTATION_LOAD)
surv = survfit(Surv(time,status)~ITHvsMutationLoad,data = survdat)
diff = survdiff(Surv(time,status)~ITHvsMutationLoad,data = survdat)


p6 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 3),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_discrete( name = "ITH(#clones), Mutation Load") +  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p6)

pdf("./Fig1E.pdf")
print(p6)
dev.off()

cphITHML = coxph(Surv(time,status) ~ clones + mutation_load, data = survdat)


survdat$ITHvsCNVLoad = strata(survdat$ITH,survdat$CNV_LOAD)
surv = survfit(Surv(time,status)~ITHvsCNVLoad,data = survdat)
diff = survdiff(Surv(time,status)~ITHvsCNVLoad,data = survdat)


p66 = ggsurv(surv) + labs(subtitle = paste("Log rank test p-value: ",format(1 - pchisq(diff$chisq,length(diff$n) - 1),digits = 3),sep = "")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + guides(linetype = F) + scale_color_manual( name = "ITH(#clones), CNV Load", values = c("#7CAE00","#F8766D", "#C77CFF","#00BFC4"), labels = c("#7CAE00" = "low, low", "#F8766D" = "low, high", "#C77CFF" = "high, low", "#00BFC4" = "high, high")) +  scale_y_continuous(labels = scales::percent, limits = c(0,1)) +  scale_x_continuous(limits = c(0,10000)) + xlab("Time (Days)") + ylab("Survival(%)")
print(p66)

pdf("./Fig1F.pdf")
print(p66)
dev.off()

cphITHCNV = coxph(Surv(time,status) ~ clones + cnv_load, data = survdat)



pg <- ggplot(data = data.frame(ITH = as.character(survdat$clones), CYT = log(survdat$CYT+1)), aes(x = ITH, y = CYT, fill = ITH)) + geom_dotplot(binaxis='y', stackdir='center',dotsize = 2.5, binwidth = 0.1,) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + scale_fill_manual(name = "ITH", values = cols[order(as.numeric(vals))], labels = vals[order(as.numeric(vals))]) + xlab("ITH") + ylab("CYT") + stat_compare_means(label.y = 10) + scale_y_continuous(breaks = seq(0,10,2))


print(pg)
ggsave(plot=pg, filename = "./Fig1G.pdf", useDingbats= F)

cmpr = list(c("low","high"))
ph <- ggplot(data = data.frame(ITH = survdat$ITH, CYT = log(survdat$CYT+1)), aes(x = ITH, y = CYT, fill = ITH)) + geom_dotplot(binaxis='y', stackdir='center',dotsize = 2.5, binwidth = 0.1) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5) + 
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), text = element_text(size = 20)) + scale_fill_manual(name = "ITH", values = c("red","blue"), labels = c("high","low")) + xlab("ITH") + ylab("CYT") + stat_compare_means(label.y = 10) +  scale_y_continuous(breaks = seq(0,10,2))


print(ph)
ggsave(plot=ph, filename = "./Fig1H.pdf", useDingbats= F)
