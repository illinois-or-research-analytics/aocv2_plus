rm(list=ls())
library(data.table); library(RPostgres)

orig_el <- fread('citing_cited_network.integer.tsv')

con <- dbConnect(RPostgres::Postgres())
rw <- dbGetQuery(con,"SELECT record_id,lower(original_paper_doi) as rw_doi FROM cleaned_rw")
orig_nl <- dbGetQuery(con,"SELECT integer_id,doi FROM exosome_dimensions_complete_nodelist")
dbDisconnect(con)

setDT(orig_nl)
setkey(orig_nl,doi)

setDT(rw)
setkey(rw,rw_doi)

remove_nodes <- orig_nl[doi %in% rw$rw_doi]

t1 <- orig_el[!(V1 %in% remove_nodes$integer_id)]
t2 <- t1[!(V2 %in% remove_nodes$integer_id)]

print('Number of edges removed after eliminating 8322 retracted publications:')
print(dim(orig_el)[1] - dim(t2)[1])

t3 <- t2[, .N, by = 'V1'][N >= 250]

t4 <- t2[!(V1 %in% t3$V1)]
t5 <- t4[!(V2 %in% t3$V1)]

print(paste0('Edges in original graph: ', dim(orig_el)[1]))
print(paste0('Edges in retraction-depleted, jc250_corrected graph: ', dim(t5)[1]))
print('Writing to file: exosome_dimensions_wedell_retraction-depleted_jc250-corrected.csv')
fwrite(t5, file='exosome_dimensions_wedell_retraction-depleted_jc250-corrected.csv')











