# This script generates configuration null models for
# citation network. I am re-using some code I wrote
# for Bradley et al (2019) that was then re-written in
# Python by the rest of the team. 

# This script was run on Jul 1, 2022 and generated 10
# shuffled networks. shuffled_ikc[1-10].csv. These
# then clustered with IKC_k10 to generate
# shuffled_ikc[1-10].clustering

# create df for citation shuffling

# wrapper script that enables passing command line args to the function perm
# George Chacko 1/5/2018
# For Scopus data I added an extra step to coerce integers to character so that R doesn't get confused. 7/4/2019

# This script takes three parameters: 
# (i) input_file: self explanatory but should be in .csv format and be a copy of a datasetxxxx table 
# from the public schema in ERNIE
# (ii) output_name_string: For example, bl_analysis_permuted_, that is used to label output files. This string
# should include _permute_ so that the Python script that calculates z_scores will work on it.
# (iii) n_permute: an integer specifying how many permutations should be performed. Typically 100-1000.
# Thus, "$ nohup Rscript permute_script.R <input_file> <output_name_string> <n_permute> &

# Command line parameters
args <- commandArgs(TRUE)
input_file <- args[1]
output_name_string <- args[2]
n_permute <- args[3]

# The permute function
perm <- function(df,output_name_string,n_permute){
# function for permuting references in Uzzi-type analysis

library(data.table); library(dplyr)
# import source file 
sorted <- fread(df,colClasses=rep('character',3))
os <- output_name_string

# loop to create background models
for (i in 1:n_permute) {

print(paste('Starting loop number',i,sep=''))
sorted[,shuffled := sample(cited),by=cited_year]

# suppress duplicates
sorted <- unique(sorted)

fwrite(sorted,file=paste(os,i,'.csv',sep=''),row.names=FALSE)
print(paste('Ended loop number',i,sep=''))
		   }
return()
}

# Call function with commmand line arguments
perm(input_file,output_name_string,n_permute);

i=0
cs1 <- fread(paste0('shuffled',i+1,'.csv'))
cs2 <- fread(paste0('shuffled',i+2,'.csv'))
cs3 <- fread(paste0('shuffled',i+3,'.csv'))
cs4 <- fread(paste0('shuffled',i+4,'.csv'))
cs5 <- fread(paste0('shuffled',i+5,'.csv'))
cs6 <- fread(paste0('shuffled',i+6,'.csv'))
cs7 <- fread(paste0('shuffled',i+7,'.csv'))
cs8 <- fread(paste0('shuffled',i+8,'.csv'))
cs9 <- fread(paste0('shuffled',i+9,'.csv'))
cs10 <- fread(paste0('shuffled',i+10,'.csv'))


