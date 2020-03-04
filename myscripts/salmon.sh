#!/bin/bash

#this directory is using bash code

cd /data/project_data/RS_RNASeq/fastq/cleanreads/


for file in BRU*H*.cl.fq

do

	echo "starting sample ${file}"
	salmon quant -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_cds_index -l A -r ${file} --validateMappings --seqBias -o /data/project_data/RS_RNASeq/salmon/allmapping/${file}

# all of our data are going to the same folder  

done


