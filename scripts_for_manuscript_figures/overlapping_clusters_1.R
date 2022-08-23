# previously named up as fig3.R
# Aug 23, 2022


# Script to compute intersections of AOC cluster
# Process to generate Fig 2 and Table 2
# George Chacko 6/30/2022

# clean workspace
rm(list=ls())
library(data.table)
library(xtable)
library(ggplot2)

# read in data from base IKC clustering and AOC clusterings

# ikc50
ikc50_k <- fread('./experiment_51/equil_IKC_50.clustering_k')
ikc50_m <- fread('./experiment_51/equil_IKC_50.clustering_mcd')
ikc50_o	<- fread('/shared/aj_manuscript_data/experiment_0/IKC_50_realignment.clustering')

# ikc40
ikc40_k <- fread('./experiment_52/equil_IKC_40.clustering_k')
ikc40_m <- fread('./experiment_52/equil_IKC_40.clustering_mcd')
ikc40_o	<- fread('/shared/aj_manuscript_data/experiment_0/IKC_40_realignment.clustering')

# ikc30
ikc30_k <- fread('./experiment_53/equil_IKC_30.clustering_k')
ikc30_m <- fread('./experiment_53/equil_IKC_30.clustering_mcd')
ikc30_o	<- fread('/shared/aj_manuscript_data/experiment_0/IKC_30_realignment.clustering')

# ikc20
ikc20_k <- fread('./experiment_54/equil_IKC_20.clustering_k')
ikc20_m <- fread('./experiment_54/equil_IKC_20.clustering_mcd')
ikc20_o	<- fread('/shared/aj_manuscript_data/experiment_0/IKC_20_realignment.clustering')

# ikc10
ikc10_k <- fread('./experiment_55/equil_IKC_10.clustering_k')
ikc10_m <- fread('./experiment_55/equil_IKC_10.clustering_mcd')
ikc10_o	<- fread('/shared/aj_manuscript_data/experiment_0/IKC_10_realignment.clustering')

#get object names
df_vec <- ls()

# import data into a list
ikc_list <- list()

for (i in 1:length(df_vec)){
ikc_list[[i]] <- get(noquote(df_vec[i]))
}
names(ikc_list) <- df_vec

# read in degree counts data 
deg <- fread('e14mrdj250cn_degree_counts.csv')

## create list of clusterings merged with degrees
fig3_list <- lapply(ikc_list,FUN=function(x) merge(x,deg,by.x='V2',by.y='integer_id'))
# tag each list element for identification
for (i in 1:length(fig3_list)){
fig3_list[[i]][,tag:=names(fig3_list)[i]]
}
# flatten list into a dataframe and write to a file
df_fig3_list <- rbindlist(fig3_list)
fwrite(df_fig3_list,file='df_fig3_list.csv')

# for pairwise intersections between clusters
# ikc10_k

# generate pairwise combinations of cluster_ids
df_nc2 <- data.frame(t(combn(unique(ikc10_k$V1),2)))

# create dataframe of intersection and union counts
intersect_df <- data.frame()
for (i in 1:dim(df_nc2)[1]){
print(i)
intersect_df[i,1] <- df_nc2[i,1] 
intersect_df[i,2] <- df_nc2[i,2] 
intersect_df[i,3] <- length(ikc10_k[V1==df_nc2[i,1]]$V2)
intersect_df[i,4] <- length(ikc10_k[V1==df_nc2[i,2]]$V2)
intersect_df[i,5] <- length(intersect(ikc10_k[V1==df_nc2[i,1]]$V2,ikc10_k[V1==df_nc2[i,2]]$V2))
intersect_df[i,6] <- length(union(ikc10_k[V1==df_nc2[i,1]]$V2,ikc10_k[V1==df_nc2[i,2]]$V2))
print("***")
}
# make colnames human readable
colnames(intersect_df) <- c("A","B","sizeA","sizeB","intersection","union")
# upgrade dataframe to data.table
setDT(intersect_df)
# add Jaccard Coefficient calculation
intersect_df[,jc:=round(intersection/union,2)]
# write to file
fwrite(intersect_df,file='intersect_df_ikc10_k)
