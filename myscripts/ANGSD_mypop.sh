myrepo="/users/a/e/aehall/EcologicalGenomics"


mkdir ${myrepo}/myresults/ANGSD

output="${myrepo}/myresults/ANGSD"

mypop="XWS"

ls /data/project_data/RS_ExomeSeq/mapping/BWA/${mypop}*sorted.rm*.bam >${output}/${mypop}_bam.list

REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# Estimating GL's and allele frequencies for all sites with ANGSD

ANGSD -b ${output}/${mypop}_bam.list \
-ref ${REF} -anc ${REF} \
-out ${output}/${mypop}_allsites \
-nThreads 1 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-setMinDepth 3 \
-minInd 2 \
-setMinDepthInd 1 \
-setMaxDepthInd 17 \ #could me PCR duplicates 
-skipTriallelic 1 \ #probably an error because youre unlikely to see more mutations at one location 
-GL 1 \ #geneotype likelihoods 
-doCounts 1 \ #generate counts at each site 
-doMajorMinor 1 \ #majorminor status equivalent across all sites. major is often ancestral allele. minor alleles are rare and often derived
-doMaf 1 \ #major minor frequencies 
-doSaf 1 \
-doHWE 1 \ #test for HW equilibrium to see that anything is majorly outside 
# -SNP_pval 1e-6