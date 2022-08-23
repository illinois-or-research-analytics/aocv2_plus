## Previously named parameterized_fig3.R
# Aug 23,2022

# Script to compute intersections of AOC cluster
# Process to generate Fig 2 and Table 2
# George Chacko 6/30/2022
# Takes a dataframe name as single parameter
# drawn from the choices in df_fig3_list.csv (see fig3_parameter_choices)
# in /shared/gc

# clean workspace
rm(list=ls())
library(data.table)
library(xtable)
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

# read in clustering data
df_fig3_list <- fread('df_fig3_list.csv')

# use a generic identifier: a475ggtf
a <- args[1]
a475ggtf <- df_fig3_list[tag==a]

# generate pairwise combinations of cluster_ids
df_nc2 <- data.frame(t(combn(unique(a475ggtf$V1),2)))

# create dataframe of intersection and union counts
intersect_df <- data.frame()
for (i in 1:dim(df_nc2)[1]){
if (i %%100 == 0){
print(i)
print("***")
}
intersect_df[i,1] <- df_nc2[i,1] 
intersect_df[i,2] <- df_nc2[i,2] 
intersect_df[i,3] <- length(a475ggtf[V1==df_nc2[i,1]]$V2)
intersect_df[i,4] <- length(a475ggtf[V1==df_nc2[i,2]]$V2)
intersect_df[i,5] <- length(intersect(a475ggtf[V1==df_nc2[i,1]]$V2,a475ggtf[V1==df_nc2[i,2]]$V2))
intersect_df[i,6] <- length(union(a475ggtf[V1==df_nc2[i,1]]$V2,a475ggtf[V1==df_nc2[i,2]]$V2))
}

# make colnames human readable
colnames(intersect_df) <- c("A","B","sizeA","sizeB","intersection","union")
# upgrade dataframe to data.table
setDT(intersect_df)
# add Jaccard Coefficient calculation
intersect_df[,jc:=round(intersection/union,2)]
# write to file
fwrite(intersect_df,file=paste0('intersect_df_',a,'.csv'))


