## Set your working directory


## Import the libraries that we're likely to need in this session
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)  ### First: BiocManager::install("vsn") AND BiocManager::install("hexbin")

## Import the counts matrix
countsTable <- read.table("RS_cds2kb_countsMatrix.txt", header=TRUE, row.names=1)
head(countsTable) # top 6 genes of our 76 samples 
dim(countsTable) # number of rows and numbers of columns 66408 number of reference transcripts 76 samples 
countsTableRound <- round(countsTable) # Need to round because DESeq wants only integers because they need to be full counts (can't have a partial count. the number there means that that is the number of reads that mapped to that specific reference eg: 26 reads mapped to sample gene MA_10001337g0010)
head(countsTableRound)

## Import the samples description table - links each sample to factors of the experimental design.
# Need the colClasses otherwise imports "day" as numeric which DESeq doesn't like, coula altneratively change to d0, d5, d10
conds <- read.delim("RS_samples.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1, colClasses=c('factor', 'factor', 'factor', 'factor')) # yes it has headers, strings are factor, the first column is row names, each of the four cols are factors because Day is listed as integers but R would have preferred them as strings. Best to use all strings in the future rather than sometihng that could be confused as a factor. 
head(conds) 
dim(conds)


## Let's see how many reads we have from each sample:
colSums(countsTableRound) #sums reads across all columns. were on the order of about 3 million reads for each sample. 
mean(colSums(countsTableRound)) #average number of reads mapped per sample
barplot(colSums(countsTableRound), las=3, cex.names=0.5,names.arg = substring(colnames(countsTableRound),1,13)) # gives an idea of how sequencing went across 76 samples. Some rae fairly low that could be of concern. Could do a fancier plot and collect by factor (eg: Hot cold)
abline(h=mean(colSums(countsTableRound)), col="blue", lwd =2)

#whats the everage number of counts per gene
rowSums(countsTableRound) # for each gene about how many counts is it
mean(rowSums(countsTableRound))
median(rowSums(countsTableRound)) #because some have really high expression and most have low expression. Mean is driving this way up. This shows dispersion across genes- differences in magnitude of expression

# what's the average number of counts per gene per sample
apply(countsTableRound,2,mean)

## Create a DESeq object and define the experimental design here with the tilde
# conds is little matrix of factors
# doesnt work from here on 
dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData = conds, design = ~ pop + day + treatment) #later we can change this to climate + data + treatment 
dim(dds)
dds  <- DESeq(dds)
resultsNames(dds)
#compares everything to the first level listed in the matrix. manually set up contrasts that you're interested in 


# Filter out genes with few reads



## Run the DESeq model to test for differential gene expression: 1) estimate size factors (per sample), 2) estimate dispersion (per gene), 3) run negative binomial glm



# List the results you've generated



# Order and list and summarize results from specific contrasts
# Here you set your adjusted p-value cutoff, can make summary tables of the number of genes differentially expressed (up- or down-regulated) for each contrast
res <- results(dds, alpha = 0.05)
res <- res[order(res$padj),]
head(res)
#automatically uses last contrast
# compares hot vs control treatments because it was the last in the model. You can change the order you put them in. In this example everything is compared to the control. log fold change is up then its high and more highly expressed in hot conditions. Makes a results matrix for each of the contrasts. 3 are higher in control and lower in the hot. 23887 are tested in total, and quite a bit had relatively low counts. In the future you could say use the who data set but discard some of the ones with low counts and improve overall yield. 
res_treatCD <- results(dds, name="treatment_D_vs_C", alpha=0.05)
res <- res[order(res_treatCD$padj),] # reorder matrix based on adjusted pvalue and call it the same thing 
head(res_treatCD)
summary(res_treatCD)
# a lot more genes are showing differential expression in this about 3% of the genes were upregulated in drought control (more highly expressed in drought)
# 434 are downregulated in drout vs. control (more highly expressed in control than drought) 
# add together upregulated and down regulated to find differential gene expression
# ~ climate + day + treatment for all then look at the specific contrasts 

##### Data visualization #####
# MA plot

plotMA(res_treatCD,ylim=c(-3,3))
#compare anything that is a postivie logfold  change will be more highly expressed (upper quad more highly expressed in hot and dry) and lower left quad more highly expressed in control. signiificant are in red and have the greatest log fold change. median counts of 24 makes sense given the distributio. overal we have relatively few reads per gene and overall we are picking up genes that are differentially expressed 
# we are relatively low on the counts. do we want to sequence more to get a stronger signature? many read means are in the hundreds. But this could be more typical of doing 3 prime RNA seq
# 
# PCA
vsd <- vst(dds, blind=FALSE) # you can force this function to oly use a certain number of genes 
data <- plotPCA(vsd, intgroup=c("climate", "treatment", "day"),returnData=TRUE)
data$treatment <- factor(data$treatment, levels=c("C","H","D"), labels = c("C","H","D"))
data$day <- factor(data$day, levels=c("0","5","10"), labels = c("0","5","10"))
ggplot(data, aes(PC1, PC2, color=climate, shape=treatment)) +
  geom_point(size=4, alpha=0.85) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  theme_minimal()

# Counts of specific top gene! (important validatition that the normalization, model is working)
d <-plotCounts(dds, gene="MA_10257300g0010", intgroup = (c("treatment","climate","day")), returnData=TRUE)
d


p <-ggplot(d, aes(x=treatment, y=count, color=day, shape=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p <- p + geom_point(position=position_jitter(w=0.3,h=0), size=3) +
  scale_x_discrete(limits=c("C","H","D"))
p

p <-ggplot(d, aes(x=treatment, y=count, shape=climate)) + 
  theme_minimal() + theme(text = element_text(size=20), panel.grid.major=element_line(colour="grey"))
p
# we told it to pull out treatment and climate and that's what will be plotted. 
#good reality check to see all these data 
# x axis is climate, y is count. its most highly expressed in the 
# the patterns look real and significant 

# Heatmap of top 20 genes sorted by pvalue
library(pheatmap)
topgenes <- head(rownames(res_treatCD),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("treatment","climate")])
pheatmap(mat, annotation_col=df)

library(pheatmap)
topgenes <- head(rownames(res_treatCD),20)
mat <- assay(vsd)[topgenes,]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[,c("treatment","climate","day")])
pheatmap(mat, annotation_col=df)
# columns are different samples and the rows are different genes. Clustered on genes and samples 
# in general climate isnt clustering, but treatment D seems to be clustered with higher expression. These top 20 genes seem to have greater expresion in the dry treatment compared to the control. Some are upregulating in response to treatment and some may be downregulating ours appear to show an upregulated response to drought ( most of the pink treatment is red. )

# in drout regulated to control genes are upregulated 

#########subsetting 
############ Try with only Day 10 data



#grep("10", names(countsTableRound), value = TRUE)
#day10countstable <- subset(countsTableRound, grep("10", names(countsTableRound), value = TRUE)) #doesn't work has to be logical



day10countstable <- countsTableRound %>% select(contains("10"))
dim(day10countstable)

dds <- DESeqDataSetFromMatrix(countData = day10countstable, colData = conds10,
                              design = ~ climate + treatment + climate:treatment)


conds10<- subset(conds, day=="10")
dim(conds10) # to separate 

# 30 samples of day 10 

