#!/bin/bash

#this directory is using bash code

cd ~/EcologicalGenomics/myresults/

mkdir fastqcTrimmed #I am creating a new directory

for file in data /data/project_data/RS_ExomeSeq/fastq/edge_fastq/XWS*fastq.gz

do

fastqc ${file} -o fastqcTrimmed/ 

#the lower case o means dont spit it to this directory put it in the folder i tell you to be in. only need to give the portion of the path that is beyond where you are. 

done


