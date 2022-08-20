# script to generate cluster_count vs tier_count
# from Akhil's data and prototype plots
# George Chacko 9/8/202
# v2 7/11/2022
# modified to correctly merge tiervclustercount_k/m data
# Aug 20, 2022
# copied to /aocv2_plus/scripts snd renamed with bl_prefix
# checked for correctness with new Baqiao Liu data

rm(list=ls())
library(data.table)
library(ggplot2)

# degree counts data
dc <- fread('e14mrdj250cn_degree_counts.csv')

# AOC_k plot
aoc_k <- fread('experiment_55/equil_IKC_10.clustering_k')
tier_k <- fread('../aj_manuscript_data/tiervclustercount_k.csv')

aoc_k_tiers <- merge(aoc_k,tier_k,by.x='V2',by.y='node_id')
colnames(aoc_k_tiers) <- c('node_id','aoc_cluster_id','V1.y',
'ikc_cluster_id','cluster_count','tier_1_count','indegree')
aoc_k_tiers[,V1.y:=NULL]
aoc_k_tiers[,indegree:=NULL]

k_merged <- merge(aoc_k_tiers,dc,by.x='node_id',by.y='integer_id')

k_merged[total_degree/1000 < 0.1,gp:='1']
k_merged[total_degree/1000 >= 0.1 & total_degree/1000 < 1,gp:='2']
k_merged[total_degree/1000 >= 1 & total_degree/1000 < 10, gp:='3']
k_merged[total_degree/1000 >= 10 & total_degree/1000 < 100, gp:='4']
k_merged[total_degree/1000 >= 100, gp:='5']

# AOC_m plot
aoc_m <- fread('experiment_55/equil_IKC_10.clustering_mcd')
tier_m <- fread('../aj_manuscript_data/tiervclustercount_mcd.csv')

aoc_m_tiers <- merge(aoc_m,tier_m,by.x='V2',by.y='node_id')
colnames(aoc_m_tiers) <- c('node_id','aoc_cluster_id','V1.y',
'ikc_cluster_id','cluster_count','tier_1_count','indegree')
aoc_m_tiers[,V1.y:=NULL]
aoc_m_tiers[,indegree:=NULL]

m_merged <- merge(aoc_m_tiers,dc,by.x='node_id',by.y='integer_id')

m_merged[total_degree/1000 < 0.1,gp:='1']
m_merged[total_degree/1000 >= 0.1 & total_degree/1000 < 1,gp:='2']
m_merged[total_degree/1000 >= 1 & total_degree/1000 < 10, gp:='3']
m_merged[total_degree/1000 >= 10 & total_degree/1000 < 100, gp:='4']
m_merged[total_degree/1000 >= 100, gp:='5']

# add tags
m_merged[,gp2:='aoc_m']
k_merged[,gp2:='aoc_k']
all_merged <- rbind(m_merged,k_merged)
fwrite(all_merged,file='aoc_all_merged_fig3.csv')
all_merged$gp2 <- factor(all_merged$gp2, levels=c('aoc_m','aoc_k'))

# Remove duplicates
t <- all_merged[,.(node_id,ikc_cluster_id,cluster_count,tier_1_count,gp,gp2)]
t1 <- unique(t[gp2=='aoc_m'])
t2 <- unique(t[gp2=='aoc_k'])
all_merged_2 <- rbind(t1,t2)

p1 <- qplot(gp,cluster_count,data=all_merged_2,geom="boxplot",group=gp,facets=.~gp2,color=gp2)+ theme_bw() +
theme(axis.text=element_text(size=18),axis.title=element_text(size=18),legend.position="none")

pdf('cluster_group.pdf')
print(p1)
dev.off()

p2 <- qplot(cluster_count,tier_1_count,data=all_merged_2,facets=gp~gp2,color=gp2) + theme_bw() +
theme(axis.text=element_text(size=18),axis.title=element_text(size=18),legend.position="none") +
theme(panel.spacing = unit(1, "lines"))

pdf('tier_cluster.pdf')
print(p2)
dev.off()





