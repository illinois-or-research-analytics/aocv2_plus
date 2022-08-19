rm(list=ls())
library(data.table)

ikc10 <- fread('/shared/aj_manuscript_data/experiment_0/IKC_10_realignment.clustering')
ikc10_c <- ikc10[,.N,by='V1']

ikc10_aoc_m <- fread('/shared/bl_aoc_redo/aug-18-rerun/correct_ikc_10_mcd.cluster')
ikc10_aoc_m_c <- ikc10_aoc_m[,.N,by='V1']

ikc10_aoc_k <- fread('/shared/bl_aoc_redo/aug-18-rerun/correct_ikc_10_k.cluster')
ikc10_aoc_k_c <- ikc10_aoc_k[,.N,by='V1']
