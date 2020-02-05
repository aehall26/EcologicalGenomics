#!/bin/bash
# We'll use this as a wrapper to run our different mapping scripts 
myrepo="/users/a/e/aehall/EcologicalGenomics"
# My population: 
mypop="XWS"
# Directory to our cleaned and paired reads: 
input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"
# define input as JUST cleaned paired reads that have my population in it 
#directory to store our outputs 
output="/data/project_data/RS_ExomeSeq/mapping"
# run mapping.sh

#to do it in a way that allows these variables to get passed into the script that's being run use the command source, then the name of the file 
source ./mapping.sh
source ./process_bam.sh
