# needs to be better documented
# George Chacko 2/18/2022

rm(list = ls())
library(data.table)
read_network <- function(file) {
	fread(file)
}

# annotate nodes by the number of references they cite
refcount <- x[, .N, by = "V1"]

# trim network progressively where t is number of references 
# ranging from 1000 to 250 as cut off

for (t in c(1000, 500, 250, 125)) {
	trims <- refcount[N >= t]
	fwrite(trims, file = paste0("citing_cited_high_references_", "t", t, ".csv"))
	trimvec <- trims[, V1]
	y <- x[!(V1 %in% trimvec)][!(V2 %in% trimvec)]
	fwrite(y, file = paste0("highref_trimmed_citing_cited_", "t", t, ".csv"))
}



