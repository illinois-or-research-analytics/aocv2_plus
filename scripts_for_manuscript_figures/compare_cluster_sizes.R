rm(list=ls())
library(data.table)
library(ggplot2)
library(xtable)
setwd('/shared/bl_aoc_redo/')

# raw clustering data are here
# /shared/bl_aoc_redo/aug-19'

# import clustering stats

# ikc50 
ikc50 <- fread('experiment_4250/k_50_correct_ikc_50_redo_k_percent_1_original_cluster_stats.csv')
ikc50_aoc_m <- fread('experiment_4250/mcd_50_correct_ikc_50_redo_k_percent_1_overlapping_cluster_stats.csv')
ikc50_aoc_k <- fread('experiment_4250/k_50_correct_ikc_50_redo_k_percent_1_overlapping_cluster_stats.csv')

# ikc40 
ikc40 <- fread('experiment_4240/k_40_correct_ikc_40_redo_k_percent_1_original_cluster_stats.csv')
ikc40_aoc_m <- fread('experiment_4240/mcd_40_correct_ikc_40_redo_k_percent_1_overlapping_cluster_stats.csv')
ikc40_aoc_k <- fread('experiment_4240/k_40_correct_ikc_40_redo_k_percent_1_overlapping_cluster_stats.csv')

# ikc30 
ikc30 <- fread('experiment_4230/k_30_correct_ikc_30_redo_k_percent_1_original_cluster_stats.csv')
ikc30_aoc_m <- fread('experiment_4230/mcd_30_correct_ikc_30_redo_k_percent_1_overlapping_cluster_stats.csv')
ikc30_aoc_k <- fread('experiment_4230/k_30_correct_ikc_30_redo_k_percent_1_overlapping_cluster_stats.csv')

# ikc20 
ikc20 <- fread('experiment_4220/k_20_correct_ikc_20_redo_k_percent_1_original_cluster_stats.csv')
ikc20_aoc_m <- fread('experiment_4220/mcd_20_correct_ikc_20_redo_k_percent_1_overlapping_cluster_stats.csv')
ikc20_aoc_k <- fread('experiment_4220/k_20_correct_ikc_20_redo_k_percent_1_overlapping_cluster_stats.csv')

# ikc10 
ikc10 <- fread('experiment_4210/k_10_correct_ikc_10_redo_k_percent_1_original_cluster_stats.csv')
ikc10_aoc_m <- fread('experiment_4210/mcd_10_correct_ikc_10_redo_k_percent_1_overlapping_cluster_stats.csv')
ikc10_aoc_k <- fread('experiment_4210/k_10_correct_ikc_10_redo_k_percent_1_overlapping_cluster_stats.csv')

t50 <- ikc50[,.(cluster_id,cluster_size)]
t50[,tag1:="o"]
t50_m <- ikc50_aoc_m[,.(cluster_id,cluster_size)]
t50_m[,tag1:="m"]
t50_k <- ikc50_aoc_k[,.(cluster_id,cluster_size)]
t50_k[,tag1:="k"]
t50_combined <- rbind(t50,t50_m,t50_k)
t50_combined[,tag2:='ikc50']

t40 <- ikc40[,.(cluster_id,cluster_size)]
t40[,tag1:="o"]
t40_m <- ikc40_aoc_m[,.(cluster_id,cluster_size)]
t40_m[,tag1:="m"]
t40_k <- ikc40_aoc_k[,.(cluster_id,cluster_size)]
t40_k[,tag1:="k"]
t40_combined <- rbind(t40,t40_m,t40_k)
t40_combined[,tag2:='ikc40']

t30 <- ikc30[,.(cluster_id,cluster_size)]
t30[,tag1:="o"]
t30_m <- ikc30_aoc_m[,.(cluster_id,cluster_size)]
t30_m[,tag1:="m"]
t30_k <- ikc30_aoc_k[,.(cluster_id,cluster_size)]
t30_k[,tag1:="k"]
t30_combined <- rbind(t30,t30_m,t30_k)
t30_combined[,tag2:='ikc30']

t20 <- ikc20[,.(cluster_id,cluster_size)]
t20[,tag1:="o"]
t20_m <- ikc20_aoc_m[,.(cluster_id,cluster_size)]
t20_m[,tag1:="m"]
t20_k <- ikc20_aoc_k[,.(cluster_id,cluster_size)]
t20_k[,tag1:="k"]
t20_combined <- rbind(t20,t20_m,t20_k)
t20_combined[,tag2:='ikc20']

t10 <- ikc10[,.(cluster_id,cluster_size)]
t10[,tag1:="o"]
t10_m <- ikc10_aoc_m[,.(cluster_id,cluster_size)]
t10_m[,tag1:="m"]
t10_k <- ikc10_aoc_k[,.(cluster_id,cluster_size)]
t10_k[,tag1:="k"]
t10_combined <- rbind(t10,t10_m,t10_k)
t10_combined[,tag2:='ikc10']

combined <- rbind(t50_combined,t40_combined,t30_combined,t20_combined,t10_combined)

ikc <- combined[tag1=='o']
aoc_m <- combined[tag1=='m']
aoc_k <- combined[tag1=='k']

ikc_m_merge <- merge(ikc,aoc_m,by.x=c('cluster_id','tag2'),by.y=c('cluster_id','tag2'))
ikc_k_merge <- merge(ikc,aoc_k,by.x=c('cluster_id','tag2'),by.y=c('cluster_id','tag2'))

p_om1 <- qplot(log(cluster_size.x),log(cluster_size.y),data=ikc_m_merge,color=tag2) + 
xlab("cluster size original") + ylab("cluster size aoc_m") + 
geom_abline(slope=1) + theme_bw()

p_om2 <- p_om1 + theme(axis.text=element_text(size=18),
      axis.title=element_text(size=20), legend.text=element_text(size=16),
      legend.title=element_blank(),legend.position=c(0.8,0.3))

pdf('bl_fig1a.pdf')
print(p_om2)
dev.off()

p_ok1 <- qplot(log(cluster_size.x),log(cluster_size.y),data=ikc_k_merge,color=tag2) +
xlab("cluster size original") + ylab("cluster size aoc_k") +
geom_abline(slope=1) + theme_bw()

p_ok2 <- p_ok1 + theme(axis.text=element_text(size=18),
      axis.title=element_text(size=20), legend.text=element_text(size=16),
      legend.title=element_blank(),legend.position=c(0.8,0.3))

pdf('bl_fig1b.pdf')
print(p_ok2)
dev.off()

# generate tables

#aoc_m
ikc_m_merge_df <- merge(ikc_m_merge[cluster_size.x==cluster_size.y][,.N,by='tag2'],
ikc_m_merge[cluster_size.x!=cluster_size.y][,.N,by='tag2'],by.x='tag2',by.y='tag2',all.x=T)
# no change to the single IKC50 cluster so edit table to correct recycling error from cbind.
ikc_m_merge_df[5,3] <- 0
colnames(ikc_m_merge_df) <- c('IKC','no_change','increase')
xtable(ikc_m_merge_df)

#aoc_k
ikc_k_merge_df <- merge(ikc_k_merge[cluster_size.x==cluster_size.y][,.N,by='tag2'],
ikc_k_merge[cluster_size.x!=cluster_size.y][,.N,by='tag2'],by.x='tag2',by.y='tag2',all.x=T)
# no change to the single IKC50 cluster so edit table to correct recycling error from cbind.
ikc_k_merge_df[5,3] <- 0
colnames(ikc_k_merge_df) <- c('IKC','no_change','increase')
xtable(ikc_k_merge_df)





