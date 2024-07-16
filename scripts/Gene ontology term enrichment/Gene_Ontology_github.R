# This script perfomrs GO enrichment on F6 data based on their bias distance from parents.

library(dplyr)
library(goseq)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(UpSetR)
library(svglite)

# First, we load in bias distance categorised data and triad (group) ids, merge them by group_id and split by genes
PCPWbiasdist_sum_categories <- read.table('PCPWbiasdist_sum_categories.tsv',header = T)

homologies <- read.csv(file="homoeologs_1_1_1_synt_and_non_synt.csv")
lnghomologies <- homologies %>%
  pivot_longer(A:D,names_to = "genome")
colnames(lnghomologies)[8] <- c("Gene")

PCPWbiasdist_homol <- PCPWbiasdist_sum_categories %>%
  inner_join(lnghomologies,by='group_id',relationship='many-to-many') %>%
  select(c(group_id,Gene,cat_over_15,cross)) %>%
  unite('category_cross',cat_over_15:cross,sep = '_')

head(PCPWbiasdist_homol)


# Load in the reference gene set for the GO enrichment analysis
all_go <- read.csv('all_GO_traes_v1-1.csv',header = T)
all_go <- all_go %>%
  inner_join(lnghomologies,by='Gene') %>%
  select(c(Gene,GO_term,group_id))
head(all_go)

# Select only genes which were considered expressed in our F6 data (mean triad TPM > 0.5)
expressed_genes <- read.table("TPMbefore_model",header = T)
expressed_genes <- unlist(unique(expressed_genes$group_id))

all_go <- subset(all_go, group_id %in% expressed_genes)
all_go <- all_go[,c('Gene','GO_term')]
dim(all_go)
length(unique(all_go$Gene)) # Number of expressed and GO annotated genes


# Create vector for gene lengths.
lengths <- read.csv("control_timecourse_gene_lengths.csv", header=T)
colnames(lengths) <- c("gene", "length")

t1 <- subset(lengths, gene %in% all_go$Gene)
head(t1)
dim(t1)

# turn into a vector called gene.lens to use with GOSeq
gene.lens <- as.numeric(t1$length)
names(gene.lens) = t1$gene
head(gene.lens)
length(gene.lens)

# Run GO enrichment analysis
GO_enriched <- data.frame(category = character(),over_represented_pvalue = numeric(),
                          under_represented_pvalue = numeric(), numDEInCat = numeric(), 
                          numInCat = numeric(), term = character(), ontology = character(),
                          over_rep_padj = numeric(),grouped_pattern=character())

assayed.genes <- as.vector(t1$gene)


for (i in (unique(PCPWbiasdist_homol$category_cross))){
  
  #now do GO stats analysis on the genes expressed in each pattern compared to all genes expressed
  #create a named binary vector for genes where one means differentially expressed and 0 means not differentially expressed
  de.genes <- (PCPWbiasdist_homol[PCPWbiasdist_homol$category_cross == i,"Gene"])
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  head(gene.vector)
  #now carry out the GOseq analysis
  pwf = nullp(gene.vector, bias.data = gene.lens, plot.fit = TRUE)
  GO.wall = goseq(pwf, gene2cat = all_go)
  # add new column with over represented GO terms padj
  GO.wall$over_rep_padj=p.adjust(GO.wall$over_represented_pvalue, method="BH")
  GO.wall[is.na(GO.wall)] <- 'N/A'
  write.table(GO.wall[GO.wall$over_rep_padj <0.01,], file = paste0("bias_dist_GO/ALL_cluster_", i, "_GOseq.tsv", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
  
  GO_enriched_cluster <- GO.wall[GO.wall$over_rep_padj <0.01 & GO.wall$ontology == "BP",]
  write.table(GO_enriched_cluster,file = paste0("bias_dist_GO/BP_cluster_", i, "_GOseq.tsv", sep = ""), sep = "\t", quote = FALSE, col.names = TRUE, row.names = F)
  head(GO_enriched_cluster)
  
  # if no enriched GO terms don't add to dataframe
  if(nrow(GO_enriched_cluster)>0) {
    
    GO_enriched_cluster$pattern <- i
    GO_enriched <- rbind(GO_enriched,GO_enriched_cluster)
  }
  #GO_enriched
  
}

#output the full table as well
write.csv(GO_enriched,'GO_enriched_full.csv')

# Check if there are overlaps between populations in DF categories

DFB_PC <- GO_enriched[GO_enriched$pattern=='DFB_PxC',c('category','term')]
DFB_PW <- GO_enriched[GO_enriched$pattern=='DFB_PxW',c('category','term')]
DFB_overlap <- DFB_PC %>%
  inner_join(DFB_PW)

DFO_a_PC <- GO_enriched[GO_enriched$pattern=='DFO_a_PxC',c('category','term')]
DFO_a_PW <- GO_enriched[GO_enriched$pattern=='DFO_a_PxW',c('category','term')]
DFO_a_overlap <- DFO_a_PC %>%
  inner_join(DFO_a_PW)

DFO_b_PC <- GO_enriched[GO_enriched$pattern=='DFO_b_PxC',c('category','term')]
DFO_b_PW <- GO_enriched[GO_enriched$pattern=='DFO_b_PxW',c('category','term')]
DFO_b_overlap <- DFO_b_PC %>%
  inner_join(DFO_b_PW)

# Check for higher overlaps by creating an upset plot out of all categories.
upsetGO <- GO_enriched[,c('category','pattern','ontology')]
upsetGO <- upsetGO %>%
  pivot_wider(names_from = 'pattern',values_from = 'ontology')
upsetGO <- as.data.frame(upsetGO)
upsetGO[is.na(upsetGO)] <- 0
upsetGO[2:11] <- lapply(upsetGO[2:11], function(x) as.numeric(gsub("BP", "1", x)))

upset(upsetGO,nsets = 10,nintersects = 100)

# Check the overlap of Conserved PxC and PxW GO terms
consPC <- GO_enriched[GO_enriched$pattern=='Conserved_PxC',]
consPW <- GO_enriched[GO_enriched$pattern=='Conserved_PxW',]
cons_overlap <- consPC %>%
  inner_join(consPW,by='category')
head(cons_overlap)



