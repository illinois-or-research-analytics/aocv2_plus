rm(list=ls())
setwd('/shared/gc/')
library(data.table); library(xtable); library(ggplot2)

singleton_k <- fread('../aj_manuscript_data/tiervclustercount_singleton_k.csv')
singleton_m <- fread('../aj_manuscript_data/tiervclustercount_singleton_mcd.csv')
ikc10 <- fread('../aj_manuscript_data/experiment_0/IKC_10_realignment.clustering')

ikc10_counts <- ikc10[,.N,by='V1'][,.(cluster_id=V1,N)]
singleton_m_counts <- singleton_m[,.N,by='cluster_id']
singleton_k_counts <- singleton_k[,.N,by='cluster_id']

t1 <- merge(ikc10_counts,singleton_m_counts,by.x='cluster_id',by.y='cluster_id')
t2 <- merge(t1,singleton_k_counts,by.x='cluster_id',by.y='cluster_id')
colnames(t2) <- c('cluster_id','ikc10','aoc_m','aoc_k')

fwrite(t2,file='singleton_one_percent_aoc.csv')

# since size increases are small and restricted only to aoc_k
t2[,perc_increase:=round(100*(aoc_k-aoc_m)/aoc_m)]

p <- qplot(perc_increase,N,data=t2[,.N,by='perc_increase']) + theme_bw() + geom_point(size=3)
p1 <- p + xlab("pct size increase") +
ylab("count of clusters") + 
theme(axis.text=element_text(size=18),axis.title=element_text(size=20))

pdf('singletons.pdf')
print(p1)
dev.off()

#t2_melt <- melt(t2,id='cluster_id')
#p <- qplot(variable,log(value),data=t2_melt,group=variable,geom="boxplot",color=variable) + 
#theme_bw() 
#p1 <- p + ylab("log cluster size") +
#theme(axis.text=element_text(size=18),axis.title=element_text(size=20),
#axis.title.x=element_blank(),legend.position='none')
#pdf('singletons.pdf')
#print(p1)
#dev.off()
#
## alternate plot
#tm2 <- melt(t2,id=c('cluster_id','ikc10'))
#
#p2 <- qplot(log(ikc10),log(value),data=tm2,facets=variable~.,color=variable) + 
#xlim(0,15) + ylim(0,15) + geom_abline(intercept=0,slope=1) + 
#xlab("log(ikc10 cluster size)") + ylab("log(AOC cluster size)") + theme_bw()
#
#p3 <- p2 + theme(axis.text=element_text(size=18),axis.title=element_text(size=20),
#legend.position='none',strip.text.y = element_text(size = 20)) 
#
#pdf('alt_singletons.pdf')
#print(p3)
#dev.off()

