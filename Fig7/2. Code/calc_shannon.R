library(entropy)
#Calculate Shannon Diversity Index


#Snyder cohort
rm(list=ls(all=TRUE))
dat<-read.table("yardena/Riaz_analysis/final_data/Clones_dat_Snyder.txt",header=T)
pts<-data.frame(unique(dat$sample_id))
pts$SDI<-0
for(a in 1:nrow(pts)){

	these_muts<-dat[dat$sample_id == pts[a,1],17]
	these_clone_counts<-table(these_muts)
	
	print(table(these_muts))
	pts[a,2]<-entropy.empirical(data.frame(these_clone_counts)[,2])
	print(pts[a,])

		     }
write.table(pts,"yardena/Riaz_analysis/final_data/Snyder_SDI_scores.txt",quote=F)

#Riaz cohort
rm(list=ls(all=TRUE))
dat<-read.table("yardena/Riaz_analysis/final_data/Clones_dat_Riaz.txt",header=T)
#Riaz data has had no previous QC, hence filter to ensure clones supported by >=2 mutations
dat<-dat[dat$samp_clust_id_count > 1, ]
pts<-data.frame(unique(dat$sample_id))
pts$SDI<-0
for(a in 1:nrow(pts)){

	these_muts<-dat[dat$sample_id == pts[a,1],5]
	these_clone_counts<-table(these_muts)
	
	print(table(these_muts))
	pts[a,2]<-entropy.empirical(data.frame(these_clone_counts)[,2])
	print(pts[a,])

		     }
write.table(pts,"yardena/Riaz_analysis/final_data/Riaz_SDI_scores.txt",quote=F)

#Hugo cohort
rm(list=ls(all=TRUE))
dat<-read.table("yardena/Riaz_analysis/final_data/Clones_dat_Hugo.txt",header=T)
pts<-data.frame(unique(dat$sample_id))
pts$SDI<-0
for(a in 1:nrow(pts)){

	these_muts<-dat[dat$sample_id == pts[a,1],4]
	these_clone_counts<-table(these_muts)
	
	print(table(these_muts))
	pts[a,2]<-entropy.empirical(data.frame(these_clone_counts)[,2])
	print(pts[a,])

		     }
write.table(pts,"yardena/Riaz_analysis/final_data/Hugo_SDI_scores.txt",quote=F)


#Van Allen cohort
rm(list=ls(all=TRUE))
dat<-read.table("yardena/Riaz_analysis/final_data/Clones_dat_VanAllen.txt",header=T)
pts<-data.frame(unique(dat$sample_id))
pts$SDI<-0
for(a in 1:nrow(pts)){

	these_muts<-dat[dat$sample_id == pts[a,1],17]
	these_clone_counts<-table(these_muts)
	
	print(table(these_muts))
	pts[a,2]<-entropy.empirical(data.frame(these_clone_counts)[,2])
	print(pts[a,])

		     }
write.table(pts,"yardena/Riaz_analysis/final_data/VanAllen_SDI_scores.txt",quote=F)

