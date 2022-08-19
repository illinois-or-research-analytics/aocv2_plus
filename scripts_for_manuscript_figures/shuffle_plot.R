rm(list=ls())
setwd('/shared/gc')
library(data.table)
library(ggplot2)

source('./aj_venv/collate_shuffled_results.R')
setwd('/shared/gc')
ikc <- fread('../aj_manuscript_data/experiment_1/k_10_totaldegree_1percent_original_cluster_stats.csv')
ikc <- ikc[,.(cluster_id,cluster_size,mcd,gp='ikc')]
shuffles <- shuffled_df[,.(cluster_id=cluster_no,cluster_size=node_count,mcd,gp=cluster_id)]
combined <- rbind(ikc,shuffles)

combined[gp=='shuffled_ikc9.clustering',gp:='s9']
combined[gp=='shuffled_ikc8.clustering',gp:='s8']
combined[gp=='shuffled_ikc7.clustering',gp:='s7']
combined[gp=='shuffled_ikc6.clustering',gp:='s6']
combined[gp=='shuffled_ikc5.clustering',gp:='s5']
combined[gp=='shuffled_ikc4.clustering',gp:='s4']
combined[gp=='shuffled_ikc3.clustering',gp:='s3']
combined[gp=='shuffled_ikc2.clustering',gp:='s2']
combined[gp=='shuffled_ikc1.clustering',gp:='s1']
combined[gp=='shuffled_ikc10.clustering',gp:='s10']

combined$gp <- factor(combined$gp, levels=c('ikc','s1','s2','s3','s4','s5',
's6','s7','s8','s9','s10'))

p <-  qplot(gp,log(cluster_size),data=combined,size=mcd) + 
scale_y_continuous(limits=c(0,15)) + theme_bw()  + 
theme(axis.text=element_text(size=18),axis.title=element_text(size=20),
legend.title=element_text(size=18))

pdf('shuffleplot.pdf')
print(p)
dev.off()








