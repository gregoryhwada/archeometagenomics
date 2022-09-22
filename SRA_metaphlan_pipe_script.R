## This set of functions generates script used to run a pipeline that takes reads from the Sequence Read Archive,
## runs them through AdapterRemoval, and then runs the trimmed reads through MetaPhlAn. The SRA_metaphlan_make_directory
## function creates a folder in the path specified with substructure for the pipeline. The functions SRA_metaphlan_pipe_paired 
## and SRA_metaphlan_pipe_single create the scripts for running pipelines for paired end reads and single end reads respectively.
## To use the pipeline, the following steps should be used.
# (1) run SRA_metaphlan_make_directory, supplying a path to a new directory for the set of reads
# (2) put the downloaded SRR_Acc_list.txt file, downloaded from the SRA run selector, in the subfolder named list
# (3) run the appropriate pipe function for either paired-end or single-end reads
# (4) pipe function will output script to run in terminal with appropriate conda environment
# The environment the script is run in will need (1) SRA Toolkit, (2) AdapterRemoval, and (3) MetaPhlAn
# Use Note: directory name should be in quotes and end with /

SRA_metaphlan_make_directory <- function(directory) {
  if(dir.exists(directory)){
    message("Directory already exisits. Aborting.")
  } else{
    dir.create(directory)
    dir.create(paste(directory,"/list",sep=""))
    dir.create(paste(directory,"/processed",sep=""))
    dir.create(paste(directory,"/metaphlan",sep=""))
  }
}


SRA_metaphlan_pipe_paired <- function(directory) {
  list <- read.csv(paste(directory,"list/","SRR_Acc_List.txt",sep=""),header = FALSE)
  n_sample <- nrow(list)
  for(i in 1:n_sample) {
      message(noquote(paste("fasterq-dump -v -S ",list[i,]," -O ",directory,sep="")))
      message(noquote(paste("adapterremoval --file1 ",directory,list[i,],"_1.fastq",
                        " --file2 ",directory,list[i,],"_2.fastq",
                        " --basename ",directory,"processed/",list[i,],
                        sep=""
      )))
      message(noquote(paste("metaphlan ",directory,"processed/",list[i,],".pair1.truncated,",
                        directory,"processed/",list[i,],".pair2.truncated",
                        " --input_type fastq --bowtie2out ",directory,"metaphlan/",list[i,],".bt2",
                        " -o ",directory,"metaphlan/",list[i,],"_profile",
                        " -t rel_ab_w_read_stats",
                        sep=""
      )))
      message("\n")
  }
}

SRA_metaphlan_pipe_single <- function(directory) {
  list <- read.csv(paste(directory,"list/","SRR_Acc_List.txt",sep=""),header = FALSE)
  n_sample <- nrow(list)
  for(i in 1:n_sample) {
    message(noquote(paste("fasterq-dump -v ",list[i,]," -O ",directory,sep="")))
    message(noquote(paste("adapterremoval --file1 ",directory,list[i,],".fastq",
                          " --basename ",directory,"processed/",list[i,],
                          sep=""
    )))
    message(noquote(paste("metaphlan ",directory,"processed/",list[i,],".truncated",
                          " --input_type fastq --bowtie2out ",directory,"metaphlan/",list[i,],".bt2",
                          " -o ",directory,"metaphlan/",list[i,],"_profile",
                          " -t rel_ab_w_read_stats",
                          sep=""
    )))
    message("\n")
  }
}

SRA_only <- function(directory) {
  list <- read.csv(paste(directory,"list/","SRR_Acc_List.txt",sep=""),header = FALSE)
  n_sample <- nrow(list)
  for(i in 1:n_sample) {
    message(noquote(paste("fasterq-dump -v ",list[i,]," -O ",directory,sep="")))
  }
}


sourcepredict_pipe <- function() {
}