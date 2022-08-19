# test script to see if R and data.table can calculate Bu's measurements
# on the exosome dataset (citing_cited_network.integer.tsv)
# The proportional calculations can be performed later
# George Chacko 2/12/2022

print(Sys.time())
print("***********")
rm(list = ls())
library(data.table)
# to prevent crashes on valhalla


x <- fread("/srv/local/shared/external/dbid/citing_cited_network.integer.tsv")
setkey(x, V2)
# setDTthreads(16)
# this became unnecessary after upgrading to the dev version of data.table
# data.table::update.dev.pkg()

# vector of all nodes that are cited at least once
bigvec <- x[, unique(V2)]
bigvec <- sort(bigvec)

# loop to generate smaller vector of nodes (avoid memory issues and speeds up)
j <-1
while (j <= length(bigvec)){
    ptm <- proc.time()
	littlevec <- bigvec[j:(j + 999)]

	# get citing articles for littlevec
	c <- x[V2 %in% littlevec] 

	# get references of littlevec
	r <- x[V1 %in% littlevec]

	# from parent dataset (x) get all edges involving neighbors of nodes in c
	cc <- x[V1 %in% c$V1 | V2 %in% c$V2]

	# from parent dataset (x) get all edges involving neighbors of nodes in r
	rr <- x[V1 %in% r$V1 | V2 %in% r$V2]
	
	# initialize df for output from nested loop below
	df <- data.frame()
	# for loop
	for (i in 1:length(littlevec)) {

		# pubs citing littlevec[i]
		citers <- cc[V2 == littlevec[i]]
		# pubs cited by littlevec[i]
		citeds <- rr[V1 == littlevec[i] ]

		# count of citing publications
		cp_level <- dim(citers)[1]

		### citing data ###		
		# For citers of a focal pub that cite each other
		citers_citers <- cc[V1 %in% citers$V1 & V2 %in% citers$V1]

		# count of links from citers of a focal pub to its other citers
		tr_citing <- dim(citers_citers)[1]

		# count of pubs citing focal pub that also cite a citer of a focal pub
		cp_r_citing_pub_nonzero <- dim(citers_citers[, .N, by = "V1"])[1]

		# count of pubs citing focal pub that do not cite a citer of a focal pub
		cp_r_citing_pub_zero <- cp_level - cp_r_citing_pub_nonzero

		### cited data ###
		# count of links from citers of a focal pub to its references 
		t1 <- cc[V1 %in% citers$V1]
		t2 <- rr[V2 %in% citeds$V2]
		cnt_link_citer_cited <- cc[V1 %in% citers$V1 & V2 %in% t2$V2]
		tr_cited <- dim(cnt_link_citer_cited)[1]
		# count of pubs citing focal pub that also cite a reference of a focal pub
		cp_r_cited_pub_nonzero <- length(unique(cnt_link_citer_cited$V1))
		# count of pubs citing focal pub that do not cite a reference of a focal pub
		cp_r_cited_pub_zero <- cp_level  - cp_r_cited_pub_nonzero

		
		tmp <- c(littlevec[i],cp_level, cp_r_citing_pub_nonzero,cp_r_citing_pub_zero, 
		    tr_citing,tr_cited, cp_r_cited_pub_nonzero, cp_r_cited_pub_zero)
		df <- rbind(df,tmp)
		colnames(df) <- c('fp', 'cp_level', 'cp_r_citing_pub_nonzero', 'cp_r_citing_pub_zero', 
		    'tr_citing','tr_cited_pub', 'cp_r_cited_pub_nonzero', 'cp_r_cited_pub_zero')
		}
j <- j+1000
# write df to csv file
fwrite(df,file=paste0('bu-6param_',j,'_',i,'csv'))
print(j)
print(i)
print("***")
print(proc.time()- ptm)
print(Sys.time())
print(" ")
}

