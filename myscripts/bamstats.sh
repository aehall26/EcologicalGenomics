#!/bin/bash
#set repo 
myrepo="/users/a/e/aehall/EcologicalGenomics"

mypop="XWS"
output="/data/project_data/RS_ExomeSeq/mapping"

echo "Num.reads R1 R2 Paired MateMapped Singletons MateMappedDiffChr" >${myrepo}/myresults/${mypop}.flagstats.txt
#name of the file mypopflagstats.txt lives in the myresults folder in my repo

for file in ${output}/BWA/${mypop}*sorted.rmdup.bam
	do 
		f=${file/.sorted.rmdup.bam/}
		name=`basename ${f}` #base name is a function and helps us rename these files 
		echo ${name} >> ${myrepo}/myresults/${mypop}.names.txt #store individual sample names
		#double arrow allows us to keep growing a single carrot would make it write it over and over again 
		samtools flagstat ${file} | awk 'NR>=6&&NR<=12 {print $1}' | column -x 
	done >>${myrepo}/myresults/${mypop}.flagstats.txt
		#pipe means take these results and pass through to next command. 
		#awk strips just the rows of data we care about. 
		#NR is number of rows greater than or equal to six and less than or equal to 12
# calculate depth of coverage from our bam files

for file in ${output}/BWA/${mypop}*sorted.rmdup.bam

	do
		samtools depth ${file} | awk '{sum+=3} END {print sum/NR}'
		#take a rolling sum of the number in the third column
	done >> ${myrepo}/myresults/${mypop}.coverage.txt
		 
