# Process to generate Fig 2 and Table 2
# modified to change field names George Chacko 7/24/2022

rm(list=ls())
library(data.table)
library(xtable)
library(ggplot2)

ikc50_k <- fread('./experiment_51/equil_IKC_50.clustering_k')
ikc50_m <- fread('./experiment_51/equil_IKC_50.clustering_mcd')
ikc50 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_50_realignment.clustering')

ikc40_k <- fread('./experiment_52/equil_IKC_40.clustering_k')
ikc40_m <- fread('./experiment_52/equil_IKC_40.clustering_mcd')
ikc40 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_40_realignment.clustering')

ikc30_k <- fread('./experiment_53/equil_IKC_30.clustering_k')
ikc30_m <- fread('./experiment_53/equil_IKC_30.clustering_mcd')
ikc30 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_30_realignment.clustering')

ikc20_k <- fread('./experiment_54/equil_IKC_20.clustering_k')
ikc20_m <- fread('./experiment_54/equil_IKC_20.clustering_mcd')
ikc20 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_20_realignment.clustering')

ikc10_k <- fread('./experiment_55/equil_IKC_10.clustering_k')
ikc10_m <- fread('./experiment_55/equil_IKC_10.clustering_mcd')
ikc10 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_10_realignment.clustering')

df_vec   <- ls()
ikc_list <- list()

for (i in 1:length(df_vec)){
ikc_list[[i]] <- get(noquote(df_vec[i]))
}
names(ikc_list) <- df_vec
ikc_list2 <- lapply(ikc_list,FUN=function(x) x[,.(no_clusters=.N),by='V2'][,.(node_count=.N),by='no_clusters'][order(no_clusters)])

for (i in 1:length(ikc_list2)){
ikc_list2[[i]][,tag:=names(ikc_list2)[i]]
}
df_ikc_list2 <- rbindlist(ikc_list2)

degree_data <- fread('e14mrdj250cn_degree_counts.csv')
df_ikc_list2[,perc:=round(100*node_count/sum(node_count),4),by='tag']
df_ikc_list2[,aoc:=paste0('aoc_',substring(tag,nchar(tag))),by=tag]


p1 <- qplot(no_clusters,log(node_count),group=tag,geom=c("point","line"),
data=df_ikc_list2[tag %like% 'ikc10'],color=tag) + theme_bw() +
xlab("clusters assigned per node") + ylab("log node count")
p2 <- p1 + theme(axis.text=element_text(size=18),axis.title=element_text(size=20),
legend.text=element_text(size=16),legend.title=element_blank(),legend.position=c(0.8,0.5))

p3 <- qplot(no_clusters,perc,group=tag,geom=c("point","line"),
data=df_ikc_list2[tag %like% 'ikc10'],color=tag) + theme_bw() +
xlab("clusters assigned per node") + ylab("percent total clustered nodes")
p4 <- p3 + theme(axis.text=element_text(size=18),axis.title=element_text(size=20),
legend.text=element_text(size=16),legend.title=element_blank(),legend.position=c(0.8,0.5))

pdf('fig2a.pdf')
print(p2)
dev.off()

pdf('fig2b.pdf')
print(p4)
dev.off()




