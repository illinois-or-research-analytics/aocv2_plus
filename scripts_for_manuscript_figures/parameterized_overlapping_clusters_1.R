## Previously named parameterized_fig3.R
# Aug 23,2022. Stripped parameter options 
# and hard coded file paths.

# Script to compute intersections of AOC clusters
# and draw new graphs with igraph

# Historical annotations
# George Chacko 6/30/2022
# Takes a dataframe name as single parameter
# drawn from the choices in df_fig3_list.csv (see fig3_parameter_choices)
# in /shared/gc

# clean workspace
rm(list=ls())
library(data.table)
library(xtable)
library(ggplot2)

# read in clustering data for ikc10_aoc_k
a <- fread('/shared/bl_aoc_redo/aug-19/correct_ikc_10_redo_k.k.cluster')
df_nc2 <- data.frame(t(combn(unique(a$V1),2)))
intersect_df <- data.frame()
for (i in 1:dim(df_nc2)[1]){
if (i %%100 == 0){
print(i)
print("***")
}
intersect_df[i,1] <- df_nc2[i,1]
intersect_df[i,2] <- df_nc2[i,2]
intersect_df[i,3] <- length(a[V1==df_nc2[i,1]]$V2)
intersect_df[i,4] <- length(a[V1==df_nc2[i,2]]$V2)
intersect_df[i,5] <- length(intersect(a[V1==df_nc2[i,1]]$V2,a[V1==df_nc2[i,2]]$V2))
intersect_df[i,6] <- length(union(a[V1==df_nc2[i,1]]$V2,a[V1==df_nc2[i,2]]$V2))
}

colnames(intersect_df) <- c("A","B","sizeA","sizeB","intersection","union")

# make colnames human readable
colnames(intersect_df) <- c("A","B","sizeA","sizeB","intersection","union")
# upgrade dataframe to data.table
setDT(intersect_df)
# add Jaccard Coefficient calculation
intersect_df[,jc:=round(intersection/union,2)]
# write to file
fwrite(intersect_df,file=paste0('intersect_df_k','ikc10_aoc_k_cluster_overlaps','.csv'))

# remove all edges where JC == 0
intersect_df2_k <- intersect_df[!jc==0]

## read in data for ikc10_aoc_m

b <- fread('/shared/bl_aoc_redo/aug-19/correct_ikc_10_redo_k.mcd.cluster')
df_nc2 <- data.frame(t(combn(unique(a$V1),2)))
intersect_df <- data.frame()
for (i in 1:dim(df_nc2)[1]){
if (i %%100 == 0){
print(i)
print("***")
}
intersect_df[i,1] <- df_nc2[i,1]
intersect_df[i,2] <- df_nc2[i,2]
intersect_df[i,3] <- length(b[V1==df_nc2[i,1]]$V2)
intersect_df[i,4] <- length(b[V1==df_nc2[i,2]]$V2)
intersect_df[i,5] <- length(intersect(b[V1==df_nc2[i,1]]$V2,b[V1==df_nc2[i,2]]$V2))
intersect_df[i,6] <- length(union(b[V1==df_nc2[i,1]]$V2,b[V1==df_nc2[i,2]]$V2))
}

colnames(intersect_df) <- c("A","B","sizeA","sizeB","intersection","union")

# make colnames human readable
colnames(intersect_df) <- c("A","B","sizeA","sizeB","intersection","union")
# upgrade dataframe to data.table
setDT(intersect_df)
# add Jaccard Coefficient calculation
intersect_df[,jc:=round(intersection/union,2)]
# write to file
fwrite(intersect_df,file=paste0('intersect_df_','ikc10_aoc_m_cluster_overlaps','.csv'))

# remove all edges where JC == 0
intersect_df2_m <- intersect_df[!jc==0]

# plotting
library(igraph)

# aoc_m graph

gm <- graph_from_data_frame(intersect_df2_m[,.(A,B,weight=jc)][weight >
 quantile(intersect_df2_m$jc,0.5)],directed=FALSE)

pdf('fig7_overlapping_m.pdf')
plot(gm,edge.color="black",vertex.color="lightblue",layout=layout_with_fr)
dev.off()

gk <- graph_from_data_frame(intersect_df2_k[,.(A,B,weight=jc)][weight >
 median(intersect_df2_k$jc)],directed=FALSE)

pdf('fig7_overlapping_k.pdf')
plot(gk,edge.color="black",vertex.color="lightblue",layout=layout_with_fr)
dev.off()
