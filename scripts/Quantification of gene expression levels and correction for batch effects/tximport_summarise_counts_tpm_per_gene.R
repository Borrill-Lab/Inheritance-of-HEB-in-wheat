# Aim is to run combine samples to gene expression level from transcript level

# Marek Glombik
#Adapted from https://mbp-tech-talks.github.io/2018-2019/intro-differential-expression/
#and https://github.com/Borrill-Lab/NAM_RNAi_Senescence/blob/main/scripts/02_tximport_summarise_counts_tpm_per_gene.R

#Steps will be:

#1: Summarise counts per gene (rather than transcript) using tximport for the studies which are to be included in the manuscript

#2: Summarise tpm per gene (rather than transcript) using tximport for the studies which are to be included in the manuscript

##### #1: Summarise counts per gene ########## 

#BiocManager::install("tximportData")
#library(tximportData)
#library(readr)
library(tximport)
library(rhdf5)

# read in pre-constructed tx2gene table (transcript to gene table)
setwd("/jic/scratch/groups/Philippa-Borrill/iwgsc_ref_seq_1.1/iwgsc_refseqv1.1_genes_2017July06/")
tx2gene <- read.csv("transcript_to_gene_refseqv1.1.csv", header=T)
head(tx2gene)


# make vector pointing to the kallisto results files   ########
setwd("/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PW/")
samples <- read.table("PWinput.txt", header=F)
samples

setwd("/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/PWkallisto/")
files <- file.path(samples$V1, "abundance.tsv", fsep ="/")
files
names(files) <- paste0(samples$V1)
head(files)
all(file.exists(files))


# read in the files and sum per gene
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
names(txi)

head(txi$abundance)

# move into directory where I will save this analysis
setwd("/jic/scratch/groups/Philippa-Borrill/Marek/eQTL-RNAseq-data-apr23/expression_analysis")

# to see counts summarised per gene
head(txi$counts)
colnames(txi$counts)

# save counts summarised per gene
write.table(txi$counts, file="PW_F6_lines_count.tsv",sep = "\t")

# to see tpm summarised per gene
#head(txi$abundance)
#colnames(txi$abundance)

# save tpm summarised per gene
write.table(txi$abundance, file="PW_F6_lines_tpm.tsv",sep = "\t")

# see lengths summarised per gene
head(txi$length)

# calculate average gene length across all samples
gene_lengths <- as.data.frame(rowMeans(txi$length))
head(gene_lengths)
colnames(gene_lengths) <- c("length")
head(gene_lengths)
#save length per gene
write.csv(gene_lengths, file="PW_F6_lines_gene_lengths.csv")

