# Script to join cleaned_rw data with exosome network to identify
# retracted articles in the network. George Chacko 3/8/2022
# needs to be run with Rscript for the readline part to work

library(RPostgres); library(data.table)

# User supplied parameters

cat("enter_host_address: ")
ha <- readLines("stdin", n=1)
cat("You entered hostname and address")
cat("\n")

cat("userid?: ")
me <- readLines("stdin", n=1)
cat("userid entered")
cat("\n")

cat("password?: ")
pw <- readLines("stdin", n=1)
cat("password entered")
cat("\n")

cat("port no?: ")
pn <- readLines("stdin", n=1)
cat("ok have port no")
cat("\n")

cat("database_name?: ")
db <- readLines("stdin", n=1)
cat("ok..will try connecting")
cat("\n")
Sys.sleep(3)

con <- dbConnect(RPostgres::Postgres(),dbname=db,host=ha, user=me, port=pn, password=pw)
# dbListTables(con)

# gets the set of nodes in the exosome network that represent publications that
# were retracted.

dbGetQuery(con,"SELECT dedcn.integer_id,dedcn.doi FROM dimensions.exosome_dimensions_complete_nodelist dedcn
INNER JOIN chackoge.cleaned_rw ccr ON lower(ccr.original_paper_doi)=lower(dedcn.doi)")

# get original paper dois joined with integer_ids
ex_nl_orig_retractions <- dbGetQuery(con,"WITH cte AS( SELECT original_paper_doi, orig_year, retraction_doi, retraction_year
FROM chackoge.cleaned_rw WHERE original_paper_doi IS NOT NULL AND original_paper_doi NOT IN ('Unavailable', 'unavailable'))
SELECT cte.original_paper_doi, cte.orig_year, cte.retraction_doi, cte.retraction_year,edcn.doi, edcn.year, edcn.integer_id
FROM cte INNER JOIN exosome_dimensions_complete_nodelist edcn ON edcn.doi=lower(cte.original_paper_doi)")

ex_nl <- dbGetQuery(con,"SELECT * FROM exosome_dimensions_complete_nodelist")

dbDisconnect(con)

setwd('~/Desktop/Retractions')
fwrite(ex_nl_orig_retractions,file='ex_orig_retractions.csv')
fwrite(ex_nl,file='ex_nl.csv')

### Alternatively
# dbSendQuery creates a remote object as far as I can tell
#ex_nl_orig_retractions <- dbSendQuery(con,
#"WITH cte AS(
#SELECT original_paper_doi, orig_year 
#FROM chackoge.cleaned_rw
#WHERE original_paper_doi IS NOT NULL
#AND original_paper_doi NOT IN ('Unavailable', 'unavailable'))
#SELECT cte.original_paper_doi, cte.orig_year, 
#	edcn.doi, edcn.year, edcn.integer_id
#FROM cte
#INNER JOIN exosome_dimensions_complete_nodelist edcn
#ON edcn.doi=lower(cte.original_paper_doi)")

# dbFetch retrieves it and should be used for chunked retrieves of large results
# enor <- dbFetch(ex_nl_orig_retractions, n= 25)

# Clean up
#dbClearResult(ex_nl_orig_retractions)
#dbDisconnect(con)

### read into R
setwd('~/Desktop/Retractions')
ex_el <- fread('citing_cited_network.integer.tsv') # exosome edgelist
ex_nl <- fread('ex_nl.csv') # exosome nodelist
ex_r <- fread('ex_orig_retractions.csv') # retraction data based on original doi matched to exosome nodelist (SQL above)
# resolve discrepancy between year(Dimensions data) and orig_year (RW) by taking the later of the two to be conservative.
ex_r[,adjusted_orig_year:=max(orig_year,year),by='integer_id']
rm(list= ls()[!(ls() %in% c('ex_r','ex_el','ex_nl'))])


# Note that 14 articles have 2 rows in ex_r on account of multiple retraction_dois (or possibly some other field)
# > ex_r[,.N,by='integer_id'][N > 1]
#  1:    4757071 2
#  2:    4806381 2
#  3:    4298595 2
#  4:   12319840 2
#  5:    1102064 2
#  6:    3161675 2
#  7:    6503172 2
#  8:    7888570 2
#  9:    6370231 2
# 10:    4701851 2
# 11:    4900446 2
# 12:   10797456 2
# 13:    8909636 2
# 14:    3960810 2

# Identify pubs that cite retracted papers (use the original paper_doi)
# %in% behaves like intersect and union in terms of not having duplicates in the output
cited_retractions <- ex_el[V2 %in% ex_r$integer_id][,.N,by='V2'][order(-N)][,.(integer_id=V2,citation_count=N)]

# For each of 4733 retracted papers that have been cited by other papers in ex_el
# Get number of citations
# Number of citations that have occurred after the retraction date
# Number of years after retraction date that citations were recorded (do we need this?)
# Whether citing papers have also been retracted
# For retractions with multiple rows take latest retraction date

dupes <- cited_retractions[integer_id %in% (ex_r[,.N,by='integer_id'][N > 1][,integer_id])]
# wc = working copy
wc <- ex_r[integer_id %in% dupes$integer_id]
wc <- ex_r[integer_id %in% dupes$integer_id] # suppress duplicate rows
wc <- wc[,.SD[which.min(retraction_year)],by=integer_id] # Alison's suggestion to use min
deduped_ex_r <- ex_r[!integer_id %in% wc$integer_id]

# use deduped ex_r for citation analysis
clean_cited_retractions <- merge(cited_retractions,deduped_ex_r,by.x='integer_id',by.y='integer_id')[order(-citation_count)]
# apply arbitrary bound of 50 citations (can be changed to some other value)
clean_cited_retractions <- clean_cited_retractions[citation_count >=50]

# a <- merge(ex_el[V2==clean_cited_retractions[i]$integer_id],ex_nl,by.x='V1',by.y='integer_id')[,.(citing=V1,retracted=V2,citing_year=year)]

# Summary data
retraction_metadata <- data.frame()
# Run a for-loop through each integer_id
for (i in 1:dim(clean_cited_retractions)[1]) {
	print(i)
	tempvec <- c(
	clean_cited_retractions[i]$integer_id,
	clean_cited_retractions[i]$citation_count,
 	clean_cited_retractions[i]$adjusted_orig_year)
 	print(tempvec)
 	retraction_metadata <- rbind(retraction_metadata,tempvec)
}
names(retraction_metadata) <- c('integer_id', 'network_citation_count', 'pub_year')

retraction_details <- data.frame()
for (i in 1:dim(clean_cited_retractions)[1]) {
a <- merge(ex_el[V2==clean_cited_retractions[i]$integer_id],ex_nl,by.x='V1',by.y='integer_id')[,.(citing=V1,retracted=V2,citing_year=year)]
a <- merge(a,clean_cited_retractions,by.x='retracted',by.y='integer_id')[,.(citing,retracted,citation_count,citing_year,retraction_year,adjusted_orig_year)]
b <- clean_cited_retractions[i]$adjusted_orig_year
a[,post_retraction_citation_period:= citing_year - b]
retraction_details <- rbind(retraction_details,a)
}

# Some citing articles are missing year of publication in ex_nl

# retraction_details[is.na(retraction_details$post_retraction_citation_period)]
#    citing retracted citing_year post_retraction_citation_period
# 1: 1672169   6150945          NA                              NA
# 2: 9462623   5815110          NA                              NA

# Fix by clunky method
retraction_details <- retraction_details[!is.na(citing_year)]

t <- unique(retraction_details[citing_year > retraction_year][,.(citation_count,.N),by='retracted'])
t[,ratio:=round(N/citation_count,2)]
library(ggplot2)

# png
png('~/repos/clustering_manuscripts/images/blog1.png')
qplot(citation_count,ratio,data=t,xlab='total_in_graph_citations',ylab='post_retraction_in_graph_citations/total_in_graph_citations',main='Post_retraction in_graph citations, exosome network \n N=236 min(total_citations=50)') + theme_bw() + geom_vline(xintercept=50,color="red") + geom_hline(yintercept=0.25,color="black")
dev.off()

#pdf
pdf('~/repos/clustering_manuscripts/images/blog1.pdf')
qplot(citation_count,ratio,data=t,xlab='total_in_graph_citations',ylab='post_retraction_in_graph_citations/total_in_graph_citations',main='Post_retraction in_graph citations, exosome network \n N=236 min(total_citations=50)') + theme_bw() + geom_vline(xintercept=50,color="red") + geom_hline(yintercept=0.25,color="black")
dev.off()



