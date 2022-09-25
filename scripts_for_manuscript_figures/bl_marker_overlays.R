# script to compare marker overlays
library(data.table)
rm(list=ls())

mkrs <- fread('/shared/gc/marker_nodes_integer_pub.csv')
ikc10 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_10_realignment.clustering')
aock10_k <- fread('/shared/bl_aoc_redo/aug-19/correct_ikc_10_redo_k.k.cluster')
aock10_m <- fread('/shared/bl_aoc_redo/aug-19/correct_ikc_10_redo_k.mcd.cluster')

# begin merging with markers and count

# ikc
# merge with markers
ikc10_m <- merge(ikc10,mkrs,by.x='V2',by.y='integer_id')
# add tag and new column names
t1 <- ikc10_m[,.(node_id=V2,cluster_id=V1,'ikc')]
# group and count
t1_counts <- t1[,.N,by='cluster_id']

# aoc_m
# merge with markers
aock10_m_m <- merge(aock10_m,mkrs,by.x='V2',by.y='integer_id')
# add tag and new column names after removing unnecessary columns
t2 <- aock10_m_m[,.(node_id=V2,aoc_m_cluster_id=V1)]
t2_counts <- t2[,.(aoc_m=.N),by='aoc_m_cluster_id']

# aoc_k
# merge with markers
aock10_k_m <- merge(aock10_k,mkrs,by.x='V2',by.y='integer_id')
# new column names after removing unnecessary columns
t3 <- aock10_k_m[,.(node_id=V2,aoc_k_cluster_id=V1)]
t3_counts <- t3[,.(aoc_k=.N),by='aoc_k_cluster_id']

temp <- merge(t2_counts,t1_counts,by.x='aoc_m_cluster_id',by.y='cluster_id',all.x=TRUE,all.y=TRUE)
t123 <- merge(t3_counts,temp,by.x='aoc_k_cluster_id',by.y='aoc_m_cluster_id',all.x=TRUE,all.y=TRUE)

dt <- t123[,.(cluster_id=aoc_k_cluster_id,ikc=N,aoc_m,aoc_k)]
# replace NA
dt[is.na(dt)] <- 0

dt[,ikc_perc:=round(100*ikc/1021)]
dt[,aoc_m_perc:=round(100*aoc_m/1021)]
dt[,aoc_k_perc:=round(100*aoc_k/1021)]

dt$aoc_m <- as.integer(dt$aoc_m)
dt$ikc_perc <- as.integer(dt$ikc_perc)
dt$aoc_m_perc <- as.integer(dt$aoc_m_perc)
dt$aoc_k_perc <- as.integer(dt$aoc_k_perc)

print(dt[aoc_k_perc > quantile(aoc_k_perc,0.9)])
library(xtable)
xtable(dt[aoc_k >0])





