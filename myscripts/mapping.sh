#!/bin/bash

# this sceipt calls the R1 and R2 reads for each individual in the population 
# bwa maps reads to the reference genome 
# options are preceeded by a - (see tutorial for specific description of the code)
ref="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"
# Write a loop to map each individual within my population 
for forward in ${input}*_R1.cl.pd.fq
#by including the wild card it says do it for every member within my population that matches my code and ends with cleaned paired .fq
do 
	reverse=${forward/_R1.cl.pd.fq/_R2.cl.pd.fq} #initialize on the forward read and then take the name of the forward read and substitute to make the name for the reverse read 
	f=${forward/_R1.cl.pd.fq/} #delete the extension 
	name=`basename ${f}` #leave us with just the name of the population. basename is a function within bash. ` means you're giving a command within what youre definin g
	bwa mem -t 1 -M ${ref} ${forward} ${reverse} > ${output}/BWA/${name}.sam
done 
