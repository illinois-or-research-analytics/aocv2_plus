rm(list=ls())
setwd('/shared/gc/aj_venv')
library(data.table)

shuffled_vec <- sort(Sys.glob("shuffled*.clustering"))

shuffled_list <- list()
for (i in 1:length(shuffled_vec)) {
shuffled_list[[i]] <- fread(shuffled_vec[i])
}
names(shuffled_list) <- shuffled_vec

shuffled_df <- data.frame()
setDT(shuffled_df)
for (i in 1:length(shuffled_list)){
temp <- shuffled_list[[i]][,.(node_count=.N), by=c('V2','V3','V4')][node_count >1]
temp[,id:=names(shuffled_list[i])]
shuffled_df <- rbind(shuffled_df,temp)
}

colnames(shuffled_df) <- c('cluster_no','mcd','modularity','node_count', 'cluster_id')
fwrite(shuffled_df,file='controlled_shuffle_summary.csv')




