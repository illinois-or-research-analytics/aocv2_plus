# revised Figure 3 in manuscript
# modified from fig2_mod.R
# George Chacko 8/20/2022

rm(list=ls())
library(data.table)
library(ggplot2)
# raw clustering from corrected aoc_m, and aoc_k scripts
setwd('/shared/bl_aoc_redo/aug-19')


## base ikc data
ikc10 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_10_realignment.clustering')
ikc20 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_20_realignment.clustering')
ikc30 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_30_realignment.clustering')
ikc40 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_40_realignment.clustering')
ikc50 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_50_realignment.clustering')

## corrected aoc_data

#aoc_m
ikc10_m <- fread('correct_ikc_10_redo_k.mcd.cluster')
ikc20_m <- fread('correct_ikc_30_redo_k.mcd.cluster')
ikc30_m <- fread('correct_ikc_30_redo_k.mcd.cluster')
ikc40_m <- fread('correct_ikc_40_redo_k.mcd.cluster')
ikc50_m <- fread('correct_ikc_50_redo_k.mcd.cluster')

#aoc_k
ikc10_k <- fread('correct_ikc_10_redo_k.k.cluster')
ikc20_k <- fread('correct_ikc_20_redo_k.k.cluster')
ikc30_k <- fread('correct_ikc_30_redo_k.k.cluster')
ikc40_k <- fread('correct_ikc_40_redo_k.k.cluster')
ikc50_k <- fread('correct_ikc_50_redo_k.k.cluster')

df_vec <- ls()
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

setwd('/shared/aocv2_plus/gc')

pdf('bl_fig2a.pdf')
print(p2)
dev.off()

pdf('bl_fig2b.pdf')
print(p4)
dev.off()




