############################################################
### Data Preparation: Metaphlan Output to Source Predict ###
############################################################

### by: Gregory Wada
## This package of R scripts transforms metaphlan taxonomic profiles into correctly formatted csv files that can be
## used as inputs for SourcePredict. The taxonomic tables can be used for additional analysis in R as well. The
## workflow benefits from having source profiles in directories based on their origin, such as from the same
## experiment on the Sequence Read Archive, or from the same environmental source.  

# attribution metaphlan - https://huttenhower.sph.harvard.edu/metaphlan/
# attribution SourcePredict - https://joss.theoj.org/papers/10.21105/joss.01540

################
## Libraries ###
################
if (!"stringi" %in% installed.packages()){
  install.packages("stringi")
}

if (!"stringr" %in% installed.packages()) {
  install.packages("stringr")
}

library("stringi") #used in writing most specific NCBI identifier
library("stringr") #used in writing most specific NCBI identifier

############################
### metaphlan_sourceprep ###
############################
## This function takes a metaphlan output created using the <-t rel_ab_w_read_stats> argument and builds data frames used to
## make two of the input files used by SourcePredict, the sample label list (sample_list) and the source table containing the
## number of reads per sample (sp_sources). The function checks if these files exist as variables in the global R environment 
## and creates them if they do not, as would be the case when running the first sample. If the files do exist, it appends the 
## information from the current sample into the running data frame. These objects will be saved in the global environment and 
## will need to be written to CSV when finished to input into SourcePredict, which can be done with write_source. Multiple 
## iterations of this function need to be run to build the data frames. The wrappers metaphlan_load, metaphlan_load_lots, or
## metaphlan_load_true automate the process. For most uses, the wrappers should be used.

# input 1: file location of metaphlan output in quotes
# input 2: sample name in quotes
# input 3: sample source in quotes (which source category does it belong to?)
# output 1: (direct output) data frame with number of reads per taxa. Row name is tax ID. Column name is sample name. 
# output 2: (global variable) running list of samples and their sources
# output 3: (global variable) running taxa table for the combined source samples

# NOTE: make sure inputs are in quotes! Alternatively, you can create the inputs as variables and update them with each run. 
# NOTE: metaphlan output needs <-t rel_ab_w_read_stats> argument. 
# Remaking metaphlan output tables is quick using the .bt2 file generated from the initial run. 
# NOTE: the default output of the function is the reduced taxa table for the given input file. Saving this to a variable can be
# used to get a reduced table with just taxa ID and count for a given sample. The sample_list and sp_sources tables will grow 
# with each iterative run and it is not necessary to save the default output anywhere if only needing these products. 
# NOTE: if a sample is entered incorrectly, it will be so in the global variable output files and can be amended through 
# standard data manipulation techniques. 

metaphlan_sourceprep <- function(mp_out, sample_name, sample_source) {
  # (1) create or append sample list in the global environment  
  if(exists("sample_list")) {
    print("Appending sample to sample_list")
    new_entry <- cbind(sample_name, sample_source)
    sample_list <<- rbind(sample_list, new_entry)
  } else {
    print("Creating sample_list")
    sample_list <- data.frame()
    new_entry <- cbind(sample_name, sample_source)
    sample_list <<- rbind(sample_list, new_entry)
  }

  # (2) create or append source list in the global environment 
  tab <- as.data.frame(read.table(mp_out, sep = "\t"))
  rownames(tab) <- tab[,2]
  tab <- as.data.frame(tab[,5,drop = FALSE])
  colnames(tab) <- sample_name
  
  if(exists("sp_sources")) {
    print("Appending data to sp_sources")
    sp_sources <- merge(sp_sources,tab,by="row.names",all='TRUE')
    sp_sources[is.na(sp_sources)] <- 0
    rownames(sp_sources) <- sp_sources[,1]
    sp_sources <<- sp_sources[,-1]
  } else {
    print("Creating sp_sources")
    sp_sources <- data.frame()
    sp_sources <<- rbind(tab)
  }
  
  return(as.data.frame(tab))
}

#################################################################
### metaphlan_load,  metphlan_load_lots,  metaphlan_load_true ###
#################################################################
## These functions allow for quick loading of multiple metaphlan profiles stored in the same directory. They will only read files
## that end in "_profile". The functions metaphlan_load and metaphlan_load_lots will assign sample names given a series name as
## input, eg. SM01, SM02, SM03, etc. The function metaphlan_load will always assign a 2 digit string for the sample number so they 
## will order alphabetically. This will not work for series of over 99 samples. In this case, the function metaphlan_load_lots can 
## be used. It will assign samples numbers with leading zeros such that all samples in the series have the same number of digits. 
## For sample names matching the original file name in the specified directory, such as their SRR number, metaphlan_load_true can 
## be used. In this case, no series number is provided as input, just the directory and the environment category. The series numbers 
## can be useful in data exploration and keeping track of samples' origins. Using the NCBI identifier can be useful in reporting. 
# input1: path to directory (in quotes) where metaphlan profiles are stored
# input2: assigned name of series (in quotes); not used in metaphlan_load_true
# input3: assigned category of environment type (in quotes)
# output: iterative run of metaphlan_sourceprep

metaphlan_load <- function(directory,series,category){
  profiles <- list.files(directory,"*_profile")
  n_sample <- length(profiles)
  sample_names <- paste(series,str_pad(1:n_sample,2,pad="0"),sep="")
  for(i in 1:n_sample) {
    metaphlan_sourceprep(paste(directory,"/",profiles[i],sep=""),sample_names[i],category) 
  }
}

metaphlan_load_lots <- function(directory,series,category){
  profiles <- list.files(directory,"*_profile")
  n_sample <- length(profiles)
  digits <- floor(log10(n_sample)) + 1
  sample_names <- paste(series,str_pad(1:n_sample,digits,pad="0"),sep="")
  for(i in 1:n_sample) {
    metaphlan_sourceprep(paste(directory,"/",profiles[i],sep=""),sample_names[i],category) 
  }
}


metaphlan_load_true <- function(directory,category){
  profiles <- list.files(directory,"*_profile")
  n_sample <- length(profiles)
  digits <- floor(log10(n_sample)) + 1
  sample_names <- sub("\\_.*", "", profiles)
  for(i in 1:n_sample) {
    metaphlan_sourceprep(paste(directory,"/",profiles[i],sep=""),sample_names[i],category) 
  }
}


metaphlan_load_sink <- function(directory,series){
  profiles <- list.files(directory,"*_profile")
  n_sample <- length(profiles)
  sample_names <- paste(series,str_pad(1:n_sample,2,pad="0"),sep="")
  for(i in 1:n_sample) {
    metaphlan_sink_combine(paste(directory,"/",profiles[i],sep=""),sample_names[i]) 
  }
}

###################
### source_bake ###
###################
## This function is used to finalize taxonomy tables when finished importing data. Data are likely being imported
## from different directories, so multiple runs of metaphlan_load are expected. SourcePredict expects NCBI taxa IDs
## as as single numerical key, while metaphlan outputs taxonomy tree information separated by the | symbol, e.g.
## 2|1224|28216|80840|80864|219181|1658672 (Ottowia oral_taxon_894). 1658672 uniquely calls this species without the
## full tree path. In some cases, the taxa names end with one or multiple | symbols. This script returns only the most
## specific taxa ID and removes trailing | symbols. This finalizes the table for export, which can be done with the
## write_source function after inspection. 
# input: sp_sources (global environment) is default; a different object can be specified.
# output: updated sp_sources data frame with most specific NCBI identifier
# Note: files will not run in AdapterRemoval if function fed MetaPhlAn files with no hits.
# Troubleshoot (above): sp_sources <- sp_sources[which(sp_sources$taxa != -1),] did not resolve this problem.

source_bake <- function(source_table=sp_sources){
  sorter <- as.character(row.names(sp_sources))
  sorter_rev <- stri_reverse(sorter)
  
  while(!all(str_sub(sorter_rev,1,1)!="|")){
    sorter_rev <- ifelse(str_sub(sorter_rev,1,1)=="|",str_sub(sorter_rev,2),sorter_rev)
    print("cleaning entries with trailing |")
  }
  
  sorter_rev <- sub("\\|.*", "", sorter_rev)
  taxa <- stri_reverse(sorter_rev)
  sp_sources$taxa <- taxa
  
  sp_sources <- aggregate(sp_sources[,which(colnames(sp_sources)!="taxa")],by=list(sp_sources$taxa),FUN=sum)
  rownames(sp_sources) <- sp_sources$Group.1
  
  sp_sources <<- sp_sources[,which(colnames(sp_sources)!="Group.1")]
  
  print("Taxa counts consolidated and given most specific NCBI ID.")
  print("sp_sources updated")
}

####################
### write_source ###
####################
## This function writes the sp_sources and sample_list data frames to csv files, which can then be used in SourcePredict.
## This function should be run after source_bake, which formats the dataframe with the correct taxa ID label. 
#input1: a path for the sp_source dataframe to be written as csv
#input2: a path for the sample_list dataframe to be written as csv
#output: csv files in the specified paths

write_source <- function(outfile_taxacount, outfile_labels){
  TAXID <- rownames(sp_sources)
  sp_sources <- cbind(TAXID,sp_sources)
  write.csv(sp_sources,outfile_taxacount,row.names = FALSE, quote = FALSE)
  
  colnames(sample_list) <- c("","labels")
  write.csv(sample_list,outfile_labels,row.names = FALSE, quote = FALSE)
}

##########################
### metaphlan_sinkprep ###
##########################
## This function prepares the csv for the sink file used by SourcePredict. The sink file is usually a single sample 
## that is compared against sources to assess likely source contributions to the sink sample. As such, this function
## asks for a single file as input, makes necessary transformations, and writes a csv for the sink in a single run. 
# input1: path to the metaphlan profile for the sink sample
# input2: path to write the sink file as csv
# input3: name of the sink sample, such as a sample number or other unique identifier

metaphlan_sinkprep <- function(in_file, out_file, sample_number) {
  t <- as.data.frame(read.table(in_file, sep = "\t"))
  t <- t[,c(2,5)]
  
  sorter <- as.character(t[,1])
  sorter_rev <- stri_reverse(sorter)
  
  while(!all(str_sub(sorter_rev,1,1)!="|")){
    sorter_rev <- ifelse(str_sub(sorter_rev,1,1)=="|",str_sub(sorter_rev,2),sorter_rev)
    print("cleaning entries with trailing |")
  }
  
  sorter_rev <- sub("\\|.*", "", sorter_rev)
  taxa <- stri_reverse(sorter_rev)
  t$taxa <- taxa
  t <- t[,-1]
  
  t <- aggregate(t[,which(colnames(t)!="taxa")],by=list(t$taxa),FUN=sum)
  rownames(t) <- t$Group.1
  colnames(t) <- c("TAXID", as.character(sample_number))
  
  write.csv(t,out_file,row.names=FALSE)
  paste("file written to",out_file)
  
  return(as.data.frame(t))
}


###

metaphlan_sink_combine <- function(mp_out, sample_name) {

  # create or append sink list in the global environment 
  tab <- as.data.frame(read.table(mp_out, sep = "\t"))
  rownames(tab) <- tab[,2]
  tab <- as.data.frame(tab[,5,drop = FALSE])
  colnames(tab) <- sample_name
  
  if(exists("sp_sink")) {
    print("Appending data to sp_sink")
    sp_sink <- merge(sp_sink,tab,by="row.names",all='TRUE')
    sp_sink[is.na(sp_sink)] <- 0
    rownames(sp_sink) <- sp_sink[,1]
    sp_sink <<- sp_sink[,-1]
  } else {
    print("Creating sp_sink")
    sp_sink <- data.frame()
    sp_sink <<- rbind(tab)
  }
  
  return(as.data.frame(tab))
}


###


sink_bake <- function(source_table=sp_sink){
  sorter <- as.character(row.names(sp_sink))
  sorter_rev <- stri_reverse(sorter)
  
  while(!all(str_sub(sorter_rev,1,1)!="|")){
    sorter_rev <- ifelse(str_sub(sorter_rev,1,1)=="|",str_sub(sorter_rev,2),sorter_rev)
    print("cleaning entries with trailing |")
  }
  
  sorter_rev <- sub("\\|.*", "", sorter_rev)
  taxa <- stri_reverse(sorter_rev)
  sp_sink$taxa <- taxa
  
  sp_sink <- aggregate(sp_sink[,which(colnames(sp_sink)!="taxa")],by=list(sp_sink$taxa),FUN=sum)
  rownames(sp_sink) <- sp_sink$Group.1
  
  sp_sink <<- sp_sink[,which(colnames(sp_sink)!="Group.1")]
  
  print("Taxa counts consolidated and given most specific NCBI ID.")
  print("sp_sink updated")
}


###

write_sink <- function(outfile_taxacount){
  TAXID <- rownames(sp_sink)
  sp_sink <- cbind(TAXID,sp_sink)
  write.csv(sp_sink,outfile_taxacount,row.names = FALSE, quote = FALSE)
}