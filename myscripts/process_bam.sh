#!/bin/bash
# this is where our output sam files are goin gto get converted into binary format (bam)
# Then were going to sort the bam files, remove the PCR duplicates, and index the duplicates

#first let's convert sam to bam and then we sort

#for f in ${output}/BWA/${mypop}*.sam

#do 
	#out=${f/.sam}
	#substitute the sam because we'll change it out for a .bam
	#sambamba-0.7.1-linux-static view -S --format=bam ${f} -o ${out}.bam #this is specific to our server -S means its a sam file. then we'll format this into a bam file
	#samtools sort ${out}.bam -o ${out}.sorted.bam #output into a sorted.bam file 
#done

#now let's remove the PCR duplicates beause we dont want the pseudo repolication that those represent

for file in ${output}/BWA/${mypop}*.sorted.bam

do
	f=${file/.sorted.bam/}
	sambamba-0.7.1-linux-static markdup -r -t 1 ${file} ${f}.sorted.rmdup.bam
	#file is defined above. #marks the PCR duplicates and removes them.
	# Assumes it's a dupiicate if you have more than one read with the exact same start and end place.
done 

#now to finish, we'll index our files. index creates coordinates on the file which allow for quick look up and retrieval 

for file in ${output}/BWA/${mypop}*.sorted.rmdup.bam

do
	samtools index ${file}
done





