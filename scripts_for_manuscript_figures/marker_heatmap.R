# script to generate a comparative heatmap of marker node
# concentration. The path to input files is hard-coded for 
# location on my Mac. George Chacko  7/4/2022
# modified on 7/21/2022

rm(list=ls())
setwd('~/Desktop/aj_manuscript')
library(tidyverse)
library(ggplot2)

mkrs <- fread('marker_nodes_integer_pub.csv')
k10 <- fread('IKC_10_realignment.clustering')
k10_k <- fread('equil_IKC_10.clustering_k')
k10_m <- fread('equil_IKC_10.clustering_mcd')

ordervec <- c("r1","r2","r3","r4","r5","r6","r7","r8","r9","r10","r11","r12","r13","r14","r15","r16")

####
#IKC
k10_mkrs <- merge(k10,mkrs,by.x='V2',by.y='integer_id',all.x=TRUE)[,.(V2 ,V1,doi)]
k10_mkrs_counts <- k10_mkrs[,.(marker_count=length(doi[!is.na(doi)]),.N),by='V1']
k10_mkrs_counts[,perc:=round(100*marker_count/1021,1)]
setkey(k10_mkrs_counts,V1)

k10_mkr_mat <- matrix(k10_mkrs_counts$perc,ncol=8,byrow=T)
colnames(k10_mkr_mat) <- paste0('c',seq_along(1:8))
rownames(k10_mkr_mat) <- paste0('r',seq_along(1:16))

k10_hm <- k10_mkr_mat %>% 
	as.data.frame() %>%
	rownames_to_column("r_id") %>%
	pivot_longer(-c(r_id), names_to = "columns", values_to = "perc")

k10_hm$r_id <- factor(k10_hm$r_id, levels=rev(ordervec))
setDT(k10_hm)
k10_hm[,gp:='ikc']

# ggplot(k10_hm, aes(x=columns,y=r_id,fill=perc)) + geom_raster() + 
# scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar") + 
# scale_x_discrete(position = "top")

####
#AOC_K

k10_k_mkrs <- merge(k10_k,mkrs,by.x='V2',by.y='integer_id',all.x=TRUE)[,.(V2 ,V1,doi)]
k10_k_mkrs_counts <- k10_k_mkrs[,.(marker_count=length(doi[!is.na(doi)]),.N),by='V1']
k10_k_mkrs_counts[,perc:=round(100*marker_count/1021,1)]
setkey(k10_k_mkrs_counts,V1)

k10_k_mkr_mat <- matrix(k10_k_mkrs_counts$perc,ncol=8,byrow=T)
colnames(k10_k_mkr_mat) <- paste0('c',seq_along(1:8))
rownames(k10_k_mkr_mat) <- paste0('r',seq_along(1:16))

k10_k_hm <- k10_k_mkr_mat %>% 
	as.data.frame() %>%
	rownames_to_column("r_id") %>%
	pivot_longer(-c(r_id), names_to = "columns", values_to = "perc")

ordervec <- c("r1",  "r2",  "r3",  "r4",  "r5",  "r6",  "r7",  "r8",  "r9", "r10", "r11", "r12", "r13", "r14", "r15", "r16")

k10_k_hm$r_id <- factor(k10_k_hm$r_id, levels=rev(ordervec))
setDT(k10_k_hm)
k10_k_hm[,gp:='aoc_k']

# ggplot(k10_k_hm, aes(x=columns,y=r_id,fill=perc)) + geom_raster() + 
# scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar") + 
# scale_x_discrete(position = "top")

#AOC_M

k10_m_mkrs <- merge(k10_m,mkrs,by.x='V2',by.y='integer_id',all.x=TRUE)[,.(V2 ,V1,doi)]
k10_m_mkrs_counts <- k10_m_mkrs[,.(marker_count=length(doi[!is.na(doi)]),.N),by='V1']
k10_m_mkrs_counts[,perc:=round(100*marker_count/1021,1)]
setkey(k10_m_mkrs_counts,V1)

k10_m_mkr_mat <- matrix(k10_m_mkrs_counts$perc,ncol=8,byrow=T)
colnames(k10_m_mkr_mat) <- paste0('c',seq_along(1:8))
rownames(k10_m_mkr_mat) <- paste0('r',seq_along(1:16))

k10_m_hm <- k10_m_mkr_mat %>% 
	as.data.frame() %>%
	rownames_to_column("r_id") %>%
	pivot_longer(-c(r_id), names_to = "columns", values_to = "perc")

ordervec <- c("r1",  "r2",  "r3",  "r4",  "r5",  "r6",  "r7",  "r8",  "r9", "r10", "r11", "r12", "r13", "r14", "r15", "r16")

k10_m_hm$r_id <- factor(k10_m_hm$r_id, levels=rev(ordervec))
setDT(k10_m_hm)
k10_m_hm[,gp:='aoc_m']

# ggplot(k10_m_hm, aes(x=columns,y=r_id,fill=perc)) + geom_raster() + 
# scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar") + 
# scale_x_discrete(position = "top")

comps <- rbind(k10_hm,k10_k_hm,k10_m_hm)
comps$gp <- factor(comps$gp, levels=c('ikc','aoc_m','aoc_k'))

# set threshold of 1% on perc
comps2 <- comps
comps2[perc < 1.0, perc:=0]

# Since omps2[perc >0,unique(r_id)] only returns rows 1-5 truncate

pdf('marker_comps.pdf')
ggplot(comps2[r_id %in% c('r1','r2','r3','r4','r5')], aes(x=columns,y=r_id,fill=perc)) + geom_tile() + 
scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar") + scale_x_discrete(position = "top") + 
facet_wrap(~gp) + xlab("column_id") + ylab("row_id") + theme(axis.text=element_text(size=11),strip.text.x = element_text(size = 18))
dev.off()

pdf('marker_comps_wide.pdf',h=3,w=6)
# ggplot(comps[r_id %in% c('r1','r2','r3','r4','r5')], aes(x=columns,y=r_id,fill=perc)) + geom_tile() + 
ggplot(comps2[r_id %in% c('r1','r2','r3','r4','r5')], aes(x=columns,y=r_id,fill=perc)) + geom_tile() + 
scale_fill_gradient2(low="darkblue", high="darkgreen", guide="colorbar") + scale_x_discrete(position = "top") + 
facet_wrap(~gp) + xlab("column_id") + ylab("row_id") + theme(axis.text=element_text(size=11),strip.text.x = element_text(size = 18))
dev.off()

# create table for paper

marker_suppl_data <- merge(merge(k10_mkrs_counts,k10_m_mkrs_counts,by.x='V1',by.y='V1'),k10_k_mkrs_counts,by.x='V1',by.y='V1')
marker_suppl_data[,perc.x:=NULL]
marker_suppl_data[,perc.y:=NULL]
colnames(marker_suppl_data) <- c('cluster_number', 'ikc mcount', 'ikc size', 'aoc_m mcount', 'aoc_m size', 'aoc_k mcount', 'aoc_k size','perc')
library(janitor); library(xtable)
marker_suppl_data <- clean_names(marker_suppl_data)
print(xtable(marker_suppl_data[aoc_k_mcount >0 ]))
fwrite(marker_suppl_data,file='~/repos/aj_manuscript/marker_suppl_data.csv')



