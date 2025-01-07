# This script will try to run the expression~genotype analysis based on each homoeologue separately and then do the triad genotyping
# It will also check the overlap of these results with cis-eQTL results to see if it works better or not
library(VariantAnnotation)
library(snpsettest)
library(tidyverse)
library(ggplot2)
library(palettetown)
library(car)
library(pheatmap)
library(lme4)
library(gridExtra)
library(VennDiagram)
library(RColorBrewer)
library(ggtern)
library(randomcoloR)
library(rdist)
library(ggpattern)

PCdiff_wide <- read.table('/Users/glombik/work/vcf_retry/gatk/PCdiff_wide',header = T)
F6_PCwhole <- read.table("/Users/glombik/work/obj1_reanalysis/F6_PCwhole",header=T,row.names = NULL,sep = "\t")
colnames(F6_PCwhole)[8:10] <- c("A_rel","B_rel","D_rel")
 
 
PC_diff_F6 <- PCdiff_wide %>%
   inner_join(F6_PCwhole)

F6_PWwhole <- read.table('/Users/glombik/work/obj1_reanalysis/F6_PWwhole',header = T,sep = '\t')

library(edgeR)

TPM_g <- read.csv(file="/Users/glombik/work/obj1_reanalysis/TPM_g_all.csv")
TPM_g <- TPM_g[TPM_g$mean_tpm > 0.5,] ## only analyse genes whose mean TPM across the 3 replicates per genotype is > 0.5 

testgene <- TPM_g[,c('gene','genotype')]
testgene <- as.data.frame(testgene)
testgene <- separate(data = testgene,col = gene, into = c("group_id","genome","tpm"),sep = "_")
testgene$group_id <- as.integer(testgene$group_id)

homologies <- read.csv(file="/Users/glombik/work/obj1_reanalysis/homoeologs_1_1_1_synt_and_non_synt.csv")

lnghomologies <- homologies %>%
  pivot_longer(A:D,names_to = "genome")
colnames(lnghomologies)[8] <- c("gene.id")

filtered_exp_genes <- testgene %>%
  left_join(lnghomologies)

filtered_exp_genes <- filtered_exp_genes[,c('gene.id','genotype')]

geneset <- as.data.frame(unique(filtered_exp_genes$gene.id))
colnames(geneset) <- c('gene.id')

rm(TPM_g)
rm(testgene)

F6countsall <- read.table("/Users/glombik/work/obj1_reanalysis/F6_lines_count.tsv",header = T)

F6countsallHC <- F6countsall[!grepl("LC$", rownames(F6countsall)), ]
rm(F6countsall)
colnames(F6countsallHC) <- gsub("\\.", "_", colnames(F6countsallHC))
natural_order <- order(colnames(F6countsallHC))
# Reorder the columns based on the natural order
F6countsallHC <- F6countsallHC[, natural_order]

genotype_list <- unique(filtered_exp_genes[,c("genotype")])

parext <- F6countsallHC
total <- geneset
for (i in seq(1, ncol(parext), by = 3)) {
  # Extract the current group of 3 columns
  columns <- parext[, i:min(i + 2, ncol(parext))]
  columns$gene.id <- rownames(columns)
  #Store the colname to find matching genotype in the genotype list
  storecolname <- unique(sub("_$", "", substr(colnames(columns), 1, 5)))
  expfilt <- filtered_exp_genes[filtered_exp_genes$genotype %in% storecolname,] %>%
    left_join(columns)
  rownames(expfilt) <- expfilt$gene.id
  expfilt <- expfilt %>%
    select(!genotype)
  geneset_i <- geneset %>%
    left_join(expfilt)
  total <- total %>%
    left_join(geneset_i)
}

total[is.na(total)] <- 0


F6filteredHCcounts <- total
rownames(F6filteredHCcounts) <- F6filteredHCcounts$gene.id

F6filteredHCcountsPC <- F6filteredHCcounts%>%
  select(starts_with('PC_')) %>%
  select(!starts_with('PC_P')) %>%
  select(!starts_with('PC_C'))
F6filteredHCcountsPW <- F6filteredHCcounts%>%
  select(starts_with('PW_')) %>%
  select(!starts_with('PW_P')) %>%
  select(!starts_with('PW_W'))

y <- F6filteredHCcountsPC
group <- factor(c(rep(1:(length(colnames(y))/3),each=3)))
y <- DGEList(counts=y,group = group)
keep <- filterByExpr(y,min.count =10,group = group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method = 'TMM')
cpmF6filteredHCcountsPC <- as.data.frame(cpm(y,log = F,normalized.lib.sizes = T))

y <- F6filteredHCcountsPW
group <- factor(c(rep(1:(length(colnames(y))/3),each=3)))
y <- DGEList(counts=y,group = group)
keep <- filterByExpr(y,min.count =10,group = group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method = 'TMM')
cpmF6filteredHCcountsPW <- as.data.frame(cpm(y,log = F,normalized.lib.sizes = T))


zcpm <- scale(cpmF6filteredHCcountsPC,center = T,scale = T)
tcpm <- t(zcpm)
cpm_pca <- prcomp(tcpm)
cpm_out <- as.data.frame(cpm_pca$x)
cpm_out$group <- sapply(strsplit(as.character(row.names(tcpm)),"_"),"[[",1)
head(cpm_out)
percentage <- round(cpm_pca$sdev / sum(cpm_pca$sdev) * 100, 2)
percentage <- paste( colnames(cpm_out), "(", paste( as.character(percentage), "%", ")", sep="") )

p <- ggplot(cpm_out,aes(x=PC1,y=PC2,color=group)) +
  geom_point(size=4.5,alpha=0.5) +
  geom_text(aes(label = rownames(cpm_out))) +
  xlab(percentage[1]) +
  ylab(percentage[2])
p

zcpm <- scale(cpmF6filteredHCcountsPW,center = T,scale = T)
tcpm <- t(zcpm)
cpm_pca <- prcomp(tcpm)
cpm_out <- as.data.frame(cpm_pca$x)
cpm_out$group <- sapply(strsplit(as.character(row.names(tcpm)),"_"),"[[",1)
head(cpm_out)
percentage <- round(cpm_pca$sdev / sum(cpm_pca$sdev) * 100, 2)
percentage <- paste( colnames(cpm_out), "(", paste( as.character(percentage), "%", ")", sep="") )

q <- ggplot(cpm_out,aes(x=PC1,y=PC2,color=group)) +
  geom_point(size=4.5,alpha=0.5) +
  geom_text(aes(label = rownames(cpm_out))) +
  xlab(percentage[1]) +
  ylab(percentage[2])
q

# save the normalised data tables
write.table(cpmF6filteredHCcountsPC,'/Users/glombik/work/vcf_retry/single_homoeolog_genotype_test/cpmF6filteredHCcountsPC',
            quote = F,row.names = T,sep = '\t')
write.table(cpmF6filteredHCcountsPW,'/Users/glombik/work/vcf_retry/single_homoeolog_genotype_test/cpmF6filteredHCcountsPW',
            quote = F,row.names = T,sep = '\t')

cpmF6filteredHCcountsPC <- read.table('/Users/glombik/work/vcf_retry/single_homoeolog_genotype_test/cpmF6filteredHCcountsPC',header = T)
cpmF6filteredHCcountsPW <- read.table('/Users/glombik/work/vcf_retry/single_homoeolog_genotype_test/cpmF6filteredHCcountsPW',header = T)


#Merge them with the separate homoeologue data - the ones that are genotyped
sep_PC_diff_F6 <- PC_diff_F6[,c('genotype','group_id','A','B','D','A_rel','B_rel','D_rel','GTfilter_A','GTfilter_B','GTfilter_D')]
sep_PC_diff_F6 <- sep_PC_diff_F6 %>%
  pivot_longer(c('A','B','D'),names_to='subgenome',values_to = 'gene') %>%
  pivot_longer(c('A_rel','B_rel','D_rel'),names_to = 'subgenome_rel',values_to = 'rel_expr') %>%
  pivot_longer(c('GTfilter_A','GTfilter_B','GTfilter_D'),names_to = 'gtfilter',values_to = 'inherited_genotype')

sep_PC_diff_F6$subgenome_rel <- gsub('_rel','',sep_PC_diff_F6$subgenome_rel)
sep_PC_diff_F6$gtfilter <- gsub('GTfilter_','',sep_PC_diff_F6$gtfilter)


sep_PC_diff_F6 <- sep_PC_diff_F6[sep_PC_diff_F6$subgenome==sep_PC_diff_F6$subgenome_rel &
                                   sep_PC_diff_F6$subgenome==sep_PC_diff_F6$gtfilter,]
sep_PC_diff_F6 <- unique(sep_PC_diff_F6)
sep_PC_diff_F6 <- sep_PC_diff_F6[!grepl("ND|ND|NE", sep_PC_diff_F6$inherited_genotype), ]
sep_PC_diff_F6 <- sep_PC_diff_F6[complete.cases(sep_PC_diff_F6$inherited_genotype),]

sep_PC_diff_F6_lessHET <- sep_PC_diff_F6[!sep_PC_diff_F6$genotype %in% c("PC_39","PC_8"),]
sep_PC_genes <- unique(sep_PC_diff_F6_lessHET[,c('gene')])


cpmF6filteredHCcountsPC$gene <- rownames(cpmF6filteredHCcountsPC)

cpm_sep_PC <- sep_PC_genes %>%
  inner_join(cpmF6filteredHCcountsPC)

cpm_sep_PC_long <- cpm_sep_PC %>%
  pivot_longer(!gene,names_to = 'f6line',values_to = 'cpm_expression')

cpm_sep_PC_long$genotype <- substr(cpm_sep_PC_long$f6line,1,5)
cpm_sep_PC_long$genotype <- gsub('_$','',cpm_sep_PC_long$genotype)

cpm_sep_PC_lessHET <- cpm_sep_PC_long %>%
  inner_join(sep_PC_diff_F6_lessHET)


sep_PW_diff_F6 <- PW_diff_F6[,c('genotype','group_id','A','B','D','A_rel','B_rel','D_rel','GTfilter_A','GTfilter_B','GTfilter_D')]
sep_PW_diff_F6 <- sep_PW_diff_F6 %>%
  pivot_longer(c('A','B','D'),names_to='subgenome',values_to = 'gene') %>%
  pivot_longer(c('A_rel','B_rel','D_rel'),names_to = 'subgenome_rel',values_to = 'rel_expr') %>%
  pivot_longer(c('GTfilter_A','GTfilter_B','GTfilter_D'),names_to = 'gtfilter',values_to = 'inherited_genotype')

sep_PW_diff_F6$subgenome_rel <- gsub('_rel','',sep_PW_diff_F6$subgenome_rel)
sep_PW_diff_F6$gtfilter <- gsub('GTfilter_','',sep_PW_diff_F6$gtfilter)


sep_PW_diff_F6 <- sep_PW_diff_F6[sep_PW_diff_F6$subgenome==sep_PW_diff_F6$subgenome_rel &
                                   sep_PW_diff_F6$subgenome==sep_PW_diff_F6$gtfilter,]
sep_PW_diff_F6 <- unique(sep_PW_diff_F6)

sep_PW_diff_F6 <- sep_PW_diff_F6[!grepl("ND|ND|NE", sep_PW_diff_F6$inherited_genotype), ]
sep_PW_diff_F6 <- sep_PW_diff_F6[complete.cases(sep_PW_diff_F6$inherited_genotype),]

sep_PW_diff_F6_lessHET <- sep_PW_diff_F6[!sep_PW_diff_F6$genotype %in% c("PW_14","PW_25","PW_38","PW_48"),]

sep_PW_genes <- unique(sep_PW_diff_F6_lessHET[,c('gene')])

cpmF6filteredHCcountsPW$gene <- rownames(cpmF6filteredHCcountsPW)

cpm_sep_PW <- sep_PW_genes %>%
  inner_join(cpmF6filteredHCcountsPW)

cpm_sep_PW_long <- cpm_sep_PW %>%
  pivot_longer(!gene,names_to = 'f6line',values_to = 'cpm_expression')

cpm_sep_PW_long$genotype <- substr(cpm_sep_PW_long$f6line,1,5)
cpm_sep_PW_long$genotype <- gsub('_$','',cpm_sep_PW_long$genotype)

cpm_sep_PW_lessHET <- cpm_sep_PW_long %>%
  inner_join(sep_PW_diff_F6_lessHET)



# Use numerical genotype values instead of characters
## load in the SNP parental-differential data
diffmap_PC <- read.table("gatk/gatk_diffmap_PC.tab",header = T)
diffmap_PW <- read.table("gatk/gatk_diffmap_PW.tab",header = T)


# get rid of depth columns
diffmap_PCnoDP <- diffmap_PC %>%
  select(!matches("DP"))

diffmap_PCnoDP[diffmap_PCnoDP=="0/0"] <- as.numeric(0)
diffmap_PCnoDP[diffmap_PCnoDP=="0/1"] <- as.numeric(1)
diffmap_PCnoDP[diffmap_PCnoDP=="1/1"] <- as.numeric(2)
diffmap_PCnoDP[diffmap_PCnoDP=="./."] <- NA


# take only F5 samples (no parents)
f6diffmap_PCnoDP <- diffmap_PCnoDP %>%
  select(!matches(c("PC_C","PC_P")))


# make f5 into long format
lngdiffmapPC <- f6diffmap_PCnoDP %>% pivot_longer(!id:ALT,names_to = "sample",values_to = 'genotype')

# Haplotype PC
PCuniq_geneGTs <- lngdiffmapPC %>%
  dplyr::count(gene.id,genotype,sample)
PCuniq_gene <- lngdiffmapPC %>%
  dplyr::count(gene.id,sample)
colnames(PCuniq_geneGTs)[4] <- "nGT"
colnames(PCuniq_gene)[3] <- "ngene"
PCuniqmerge <- PCuniq_geneGTs %>%
  inner_join(PCuniq_gene)
PCuniqwide <- PCuniqmerge %>%
  pivot_wider(names_from = genotype,values_from = nGT)
colnames(PCuniqwide)[4:7] <- c('ref','het','alt','naa')

PCuniqwide[is.na(PCuniqwide)] <- 0
PCuniqwide$het.ratio <- PCuniqwide$het/PCuniqwide$ngene
PCuniqwide$ref.ratio <- PCuniqwide$ref/PCuniqwide$ngene
PCuniqwide$alt.ratio <- PCuniqwide$alt/PCuniqwide$ngene
PCuniqwide$naa.ratio <- PCuniqwide$naa/PCuniqwide$ngene


PCuniqwide$GTfilter <- "ND" # Not Distinguishable
PCuniqwide$GTfilter[PCuniqwide$ref.ratio== 1 | 
                      (PCuniqwide$ref.ratio + PCuniqwide$naa.ratio == 1 & 
                         PCuniqwide$naa.ratio > 0 & PCuniqwide$ref.ratio > 0)] <- 0
PCuniqwide$GTfilter[PCuniqwide$alt.ratio== 1 | 
                      (PCuniqwide$alt.ratio + PCuniqwide$naa.ratio == 1 & 
                         PCuniqwide$naa.ratio > 0 & PCuniqwide$alt.ratio > 0)] <- 2
PCuniqwide$GTfilter[PCuniqwide$het.ratio== 1 | 
                      (PCuniqwide$het.ratio + PCuniqwide$naa.ratio == 1 & 
                         PCuniqwide$naa.ratio > 0 & PCuniqwide$het.ratio > 0)] <- 1

#check HET ratio compared to other SNPs
hetPC <- PCuniqwide %>%
  group_by(sample,GTfilter) %>%
  dplyr::count(GTfilter) %>%
  group_by(sample) %>%
  mutate(total_count=sum(n),
         ratio_het=sum(n[GTfilter==1]) / total_count*100)
ggplot(hetPC,aes(sample,ratio_het)) +
  geom_point()
# two samples display extremely higher HET from other PxC lines (PC_39 and PC_8)

PCdiff_homol <- PCuniqwide %>%
  inner_join(lnghomologies) %>%
  dplyr::select(!(ngene:naa.ratio))

colnames(PCdiff_homol)[1:2] <- c('gene','genotype')
PCdiff_homol$genotype <- gsub('_$','',PCdiff_homol$genotype)
write.table(PCdiff_homol,file='PCdiff_homol.tsv',quote = F,row.names = F,sep = '\t')

PC_cpm_reSNP <- cpm_sep_PC_lessHET %>%
  inner_join(PCdiff_homol)
PC_cpm_reSNP <- PC_cpm_reSNP[PC_cpm_reSNP$GTfilter!='ND',]
PC_cpm_reSNP$GTfilter <- as.numeric(PC_cpm_reSNP$GTfilter)

write.table(PC_cpm_reSNP,file='PC_cpm_reSNP.tsv',quote = F,row.names = F,sep = '\t')


# Run model where at least > 80 % samples are present 
PC_cpm_reSNP_lm <- PC_cpm_reSNP %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(n()>120) %>%
  dplyr::summarise(pval = Anova(lm(formula=cpm_expression~GTfilter))[1,4])
PC_cpm_reSNP_lm$adjustP <- p.adjust(PC_cpm_reSNP_lm$pval,method = "BH")
PC_cpm_reSNP_lm$gt_sig <- T
PC_cpm_reSNP_lm <- PC_cpm_reSNP_lm[complete.cases(PC_cpm_reSNP_lm),]
PC_cpm_reSNP_lm[PC_cpm_reSNP_lm$adjustP > 0.05,]$gt_sig <- F
PC_cpm_reSNP_lm <- PC_cpm_reSNP_lm[-c(2)]
table(PC_cpm_reSNP_lm$gt_sig)
# FALSE  TRUE 
# 498  1093 

# Get the overall number for triads analysed, passed, failed
overallPC_cpm_reSNP <- PC_cpm_reSNP_lm %>%
  inner_join(PC_cpm_reSNP)
length(unlist(unique(overallPC_cpm_reSNP[,7])))
# 1306
length(unlist(unique(overallPC_cpm_reSNP[overallPC_cpm_reSNP$gt_sig==TRUE,7])))
# 980
length(unlist(unique(overallPC_cpm_reSNP[overallPC_cpm_reSNP$gt_sig==FALSE,7])))
# 461

#Numbers dont add up so there are triads where we have multiple homoeologs analysed and one of them is True, another is False
doublehitsPC <- overallPC_cpm_reSNP %>%
  dplyr::select(c('gt_sig','group_id')) %>%
  distinct() %>%
  group_by(group_id) %>%
  summarise(n_group_id=n())
table(doublehitsPC$n_group_id)

# PW
# get rid of depth columns
diffmap_PWnoDP <- diffmap_PW %>%
  select(!matches("DP"))

diffmap_PWnoDP[diffmap_PWnoDP=="0/0"] <- as.numeric(0)
diffmap_PWnoDP[diffmap_PWnoDP=="0/1"] <- as.numeric(1)
diffmap_PWnoDP[diffmap_PWnoDP=="1/1"] <- as.numeric(2)
diffmap_PWnoDP[diffmap_PWnoDP=="./."] <- NA


# take only F5 samples (no parents)
f6diffmap_PWnoDP <- diffmap_PWnoDP %>%
  select(!matches(c("PW_W","PW_P")))


# make f5 into long format
lngdiffmapPW <- f6diffmap_PWnoDP %>% pivot_longer(!id:ALT,names_to = "sample",values_to = 'genotype')

# Haplotype PW
PWuniq_geneGTs <- lngdiffmapPW %>%
  dplyr::count(gene.id,genotype,sample)
PWuniq_gene <- lngdiffmapPW %>%
  dplyr::count(gene.id,sample)
colnames(PWuniq_geneGTs)[4] <- "nGT"
colnames(PWuniq_gene)[3] <- "ngene"
PWuniqmerge <- PWuniq_geneGTs %>%
  inner_join(PWuniq_gene)
PWuniqwide <- PWuniqmerge %>%
  pivot_wider(names_from = genotype,values_from = nGT)
colnames(PWuniqwide)[4:7] <- c('ref','het','alt','naa')

PWuniqwide[is.na(PWuniqwide)] <- 0
PWuniqwide$het.ratio <- PWuniqwide$het/PWuniqwide$ngene
PWuniqwide$ref.ratio <- PWuniqwide$ref/PWuniqwide$ngene
PWuniqwide$alt.ratio <- PWuniqwide$alt/PWuniqwide$ngene
PWuniqwide$naa.ratio <- PWuniqwide$naa/PWuniqwide$ngene


PWuniqwide$GTfilter <- "ND" # Not Distinguishable
PWuniqwide$GTfilter[PWuniqwide$ref.ratio== 1 | 
                      (PWuniqwide$ref.ratio + PWuniqwide$naa.ratio == 1 & 
                         PWuniqwide$naa.ratio > 0 & PWuniqwide$ref.ratio > 0)] <- 0
PWuniqwide$GTfilter[PWuniqwide$alt.ratio== 1 | 
                      (PWuniqwide$alt.ratio + PWuniqwide$naa.ratio == 1 & 
                         PWuniqwide$naa.ratio > 0 & PWuniqwide$alt.ratio > 0)] <- 2
PWuniqwide$GTfilter[PWuniqwide$het.ratio== 1 | 
                      (PWuniqwide$het.ratio + PWuniqwide$naa.ratio == 1 & 
                         PWuniqwide$naa.ratio > 0 & PWuniqwide$het.ratio > 0)] <- 1

#check HET ratio compared to other SNPs
hetPW <- PWuniqwide %>%
  group_by(sample,GTfilter) %>%
  dplyr::count(GTfilter) %>%
  group_by(sample) %>%
  mutate(total_count=sum(n),
         ratio_het=sum(n[GTfilter==1]) / total_count*100)
ggplot(hetPW,aes(sample,ratio_het)) +
  geom_point()
# four samples display extremely high HET from othe PxW samples (PW_14, PW_25, PW_38 and PW_48)


PWdiff_homol <- PWuniqwide %>%
  inner_join(lnghomologies) %>%
  dplyr::select(!(ngene:naa.ratio))

colnames(PWdiff_homol)[1:2] <- c('gene','genotype')
PWdiff_homol$genotype <- gsub('_$','',PWdiff_homol$genotype)
write.table(PWdiff_homol,file='PWdiff_homol.tsv',quote = F,row.names = F,sep = '\t')

PW_cpm_reSNP <- cpm_sep_PW_lessHET %>%
  inner_join(PWdiff_homol)
PW_cpm_reSNP <- PW_cpm_reSNP[PW_cpm_reSNP$GTfilter!='ND',]
PW_cpm_reSNP$GTfilter <- as.numeric(PW_cpm_reSNP$GTfilter)

write.table(PW_cpm_reSNP,file='PW_cpm_reSNP.tsv',quote = F,row.names = F,sep = '\t')


# Run model where at least > 80 % samples are present 
PW_cpm_reSNP_lm <- PW_cpm_reSNP %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(n()>120) %>%
  dplyr::summarise(pval = Anova(lm(formula=cpm_expression~GTfilter))[1,4])
PW_cpm_reSNP_lm$adjustP <- p.adjust(PW_cpm_reSNP_lm$pval,method = "BH")
PW_cpm_reSNP_lm$gt_sig <- T
PW_cpm_reSNP_lm <- PW_cpm_reSNP_lm[complete.cases(PW_cpm_reSNP_lm),]
PW_cpm_reSNP_lm[PW_cpm_reSNP_lm$adjustP > 0.05,]$gt_sig <- F
PW_cpm_reSNP_lm <- PW_cpm_reSNP_lm[-c(2)]
table(PW_cpm_reSNP_lm$gt_sig)
# FALSE  TRUE 
# 219  484 

# Get the overall number for triads analysed
overallPW_cpm_reSNP <- PW_cpm_reSNP_lm %>%
  inner_join(PW_cpm_reSNP)

length(unlist(unique(overallPW_cpm_reSNP[,c('group_id')])))
# 658
length(unlist(unique(overallPW_cpm_reSNP[overallPW_cpm_reSNP$gt_sig==TRUE,7])))
# 464
length(unlist(unique(overallPW_cpm_reSNP[overallPW_cpm_reSNP$gt_sig==FALSE,7])))
# 217

#Numbers dont add up so there are triads where we have multiple homoeologs analysed and one of them is True, another is False
doublehitsPW <- overallPW_cpm_reSNP %>%
  dplyr::select(c('gt_sig','group_id')) %>%
  distinct() %>%
  group_by(group_id) %>%
  summarise(n_group_id=n())
table(doublehitsPW$n_group_id)

#Export overall significant homoeologues + triad info
short_sig_sep_PC <- unique(overallPC_cpm_reSNP[overallPC_cpm_reSNP$gt_sig==T,c('gene','group_id')])
short_sig_sep_PW <- unique(overallPW_cpm_reSNP[overallPW_cpm_reSNP$gt_sig==T,c('gene','group_id')])

write.table(short_sig_sep_PC,file='/Users/glombik/work/vcf_retry/significant_homoeologs_triadsPC.tsv',sep = '\t',quote = F,row.names = F)
write.table(short_sig_sep_PW,file='/Users/glombik/work/vcf_retry/significant_homoeologs_triadsPW.tsv',sep = '\t',quote = F,row.names = F)

# look at overlap of associated homoeologues and expr with bias dist results
PCPWgraph_sumbartab <- read.table('PCPWbiasdist_sum_categories.tsv',header = T)
PCsum <- PCPWgraph_sumbartab[PCPWgraph_sumbartab$cross=='PxC',]
PWsum <- PCPWgraph_sumbartab[PCPWgraph_sumbartab$cross=='PxW',]



join_bias_short_cpmPC <- PC_cpm_reSNP %>%
  inner_join(PC_cpm_reSNP_lm) %>%
  dplyr::select(c('group_id','gt_sig'))
join_bias_short_cpmPC <- unique(join_bias_short_cpmPC)
join_bias_short_cpmPC <- join_bias_short_cpmPC %>%
  inner_join(PCsum,relationship = 'many-to-many')
join_bias_short_cpmPC <- unique(join_bias_short_cpmPC)
join_fail_bias_short_cpmPC <- join_bias_short_cpmPC[join_bias_short_cpmPC$gt_sig==FALSE,]
join_bias_short_cpmPC <- join_bias_short_cpmPC[join_bias_short_cpmPC$gt_sig==TRUE,]
join_fail_bias_short_cpmPC <- join_fail_bias_short_cpmPC %>%
  anti_join(join_bias_short_cpmPC,by=c('group_id'))
table(join_bias_short_cpmPC$cat_over_15)
table(join_fail_bias_short_cpmPC$cat_over_15)
table(PCsum$cat_over_15)

#What about the rest? Did they not enter model?

overlap_PC_reSNP <- PC_cpm_reSNP %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(n()<120) %>%
  dplyr::select(c('group_id'))
overlap_PC_reSNP <- overlap_PC_reSNP[,2]
overlap_PC_reSNP <- unique(overlap_PC_reSNP)
overlap_PC_reSNP <- overlap_PC_reSNP %>%
  inner_join(PCsum,relationship = 'many-to-many') %>%
  anti_join(join_bias_short_cpmPC) %>%
  anti_join(join_fail_bias_short_cpmPC)
overlap_PC_reSNP <- unique(overlap_PC_reSNP)
table(overlap_PC_reSNP$cat_over_15)

# Still some missing...is it unclear genotype?
unclear_PC <- cpm_sep_PC_lessHET %>%
  inner_join(PCdiff_homol) %>%
  dplyr::select(c('group_id','gene'))
unclear_PC <- unique(unclear_PC)
unclear_PC <- unclear_PC[,1]
unclear_PC <- unclear_PC %>%
  anti_join(join_bias_short_cpmPC) %>%
  anti_join(join_fail_bias_short_cpmPC) %>%
  anti_join(overlap_PC_reSNP) %>%
  inner_join(PCsum,relationship = 'many-to-many')
unclear_PC <- unique(unclear_PC)
table(unclear_PC$cat_over_15)

# Still some missing...does the total SNP count fit?
SNPunclear_PC <- cpm_sep_PC_lessHET %>%
  inner_join(PCdiff_homol) %>%
  dplyr::select(c('group_id','gene'))
SNPunclear_PC <- unique(SNPunclear_PC)
SNPunclear_PC <- SNPunclear_PC[,1]
SNPunclear_PC <- SNPunclear_PC %>%
  inner_join(PCsum,relationship = 'many-to-many')
SNPunclear_PC <- unique(SNPunclear_PC)
table(SNPunclear_PC$cat_over_15)

#Now for PW

join_bias_short_cpmPW <- PW_cpm_reSNP %>%
  inner_join(PW_cpm_reSNP_lm) %>%
  dplyr::select(c('group_id','gt_sig'))
join_bias_short_cpmPW <- unique(join_bias_short_cpmPW)
join_bias_short_cpmPW <- join_bias_short_cpmPW %>%
  inner_join(PWsum,relationship = 'many-to-many')
join_bias_short_cpmPW <- unique(join_bias_short_cpmPW)
join_fail_bias_short_cpmPW <- join_bias_short_cpmPW[join_bias_short_cpmPW$gt_sig==FALSE,]
join_bias_short_cpmPW <- join_bias_short_cpmPW[join_bias_short_cpmPW$gt_sig==TRUE,]
join_fail_bias_short_cpmPW <- join_fail_bias_short_cpmPW %>%
  anti_join(join_bias_short_cpmPW,by=c('group_id'))
table(join_bias_short_cpmPW$cat_over_15)
table(join_fail_bias_short_cpmPW$cat_over_15)
table(PWsum$cat_over_15)

#What about the rest? Did they not enter model?

overlap_PW_reSNP <- PW_cpm_reSNP %>%
  dplyr::group_by(gene) %>%
  dplyr::filter(n()<120) %>%
  dplyr::select(c('group_id'))
overlap_PW_reSNP <- overlap_PW_reSNP[,2]
overlap_PW_reSNP <- unique(overlap_PW_reSNP)
overlap_PW_reSNP <- overlap_PW_reSNP %>%
  inner_join(PWsum,relationship = 'many-to-many') %>%
  anti_join(join_bias_short_cpmPW) %>%
  anti_join(join_fail_bias_short_cpmPW)
overlap_PW_reSNP <- unique(overlap_PW_reSNP)
table(overlap_PW_reSNP$cat_over_15)

# Still some missing...is it unclear genotype?
unclear_PW <- cpm_sep_PW_lessHET %>%
  inner_join(PWdiff_homol) %>%
  dplyr::select(c('group_id','gene'))
unclear_PW <- unique(unclear_PW)
unclear_PW <- unclear_PW[,1]
unclear_PW <- unclear_PW %>%
  anti_join(join_bias_short_cpmPW) %>%
  anti_join(join_fail_bias_short_cpmPW) %>%
  anti_join(overlap_PW_reSNP) %>%
  inner_join(PWsum,relationship = 'many-to-many')
unclear_PW <- unique(unclear_PW)
table(unclear_PW$cat_over_15)

# Still some missing...does the total SNP count fit?
SNPunclear_PW <- cpm_sep_PW_lessHET %>%
  inner_join(PWdiff_homol) %>%
  dplyr::select(c('group_id','gene'))
SNPunclear_PW <- unique(SNPunclear_PW)
SNPunclear_PW <- SNPunclear_PW[,1]
SNPunclear_PW <- SNPunclear_PW %>%
  inner_join(PWsum,relationship = 'many-to-many')
SNPunclear_PW <- unique(SNPunclear_PW)
table(SNPunclear_PW$cat_over_15)

# Now analyse whether the homoeolog associated with expression corresponds with the change in the triad direction
# Extract just the genes that were found associated + the homoeolog info and check it through Arun's tpm table

checkdirPC <- PC_cpm_reSNP_lm[PC_cpm_reSNP_lm$gt_sig==TRUE,]
checkdirPC$homoeologue <- substr(checkdirPC$gene,9,9)

genegroup <- lnghomologies
colnames(genegroup)[8] <- 'gene'

checkdirPC <- checkdirPC %>%
  left_join(genegroup)


checkdirPC <- checkdirPC[,c('group_id','homoeologue')]

checkdirPC <- checkdirPC %>%
  group_by(group_id) %>%
  summarise(associated_homoeologue = paste(homoeologue, collapse = ","))

genotPC <- unique(PC_cpm_reSNP[,c('group_id','subgenome','inherited_genotype','genotype')])
genotPC <- genotPC %>%
  arrange(group_id,subgenome) %>%
  group_by(group_id,genotype) %>%
  summarise(SNP_for_subgenome_collapsed = paste(subgenome,collapse = ','),inherited_genotype_collapsed=paste(inherited_genotype,collapse = ','))

updated_F6_PCwhole <- F6_PCwhole %>%
  inner_join(checkdirPC) %>%
  inner_join(genotPC) %>%
  dplyr::select(!cv_parent:tpm_parent)
updated_F6_PCwhole <- updated_F6_PCwhole %>%
  pivot_wider(names_from = which_parent,values_from = distance_to_parent)


# PW

checkdirPW <- PW_cpm_reSNP_lm[PW_cpm_reSNP_lm$gt_sig==TRUE,]
checkdirPW$homoeologue <- substr(checkdirPW$gene,9,9)

genegroup <- lnghomologies
colnames(genegroup)[8] <- 'gene'

checkdirPW <- checkdirPW %>%
  left_join(genegroup)


checkdirPW <- checkdirPW[,c('group_id','homoeologue')]

checkdirPW <- checkdirPW %>%
  group_by(group_id) %>%
  summarise(associated_homoeologue = paste(homoeologue, collapse = ","))

genotPW <- unique(PW_cpm_reSNP[,c('group_id','subgenome','inherited_genotype','genotype')])
genotPW <- genotPW %>%
  arrange(group_id,subgenome) %>%
  group_by(group_id,genotype) %>%
  summarise(SNP_for_subgenome_collapsed = paste(subgenome,collapse = ','),inherited_genotype_collapsed=paste(inherited_genotype,collapse = ','))

updated_F6_PWwhole <- F6_PWwhole %>%
  inner_join(checkdirPW) %>%
  inner_join(genotPW) %>%
  dplyr::select(!cv_parent:tpm_parent)

updated_F6_PWwhole <- updated_F6_PWwhole %>%
  pivot_wider(names_from = which_parent,values_from = distance_to_parent)



write.table(updated_F6_PCwhole,'/Users/glombik/work/vcf_retry/updated_F6_PCwhole',quote = F,row.names = F,sep = '\t')
write.table(updated_F6_PWwhole,'/Users/glombik/work/vcf_retry/updated_F6_PWwhole',quote = F,row.names = F,sep = '\t')

#Create supplementary tables with all the necessary information inside
PCdistcheck <- F6_PCwhole %>%
  pivot_wider(values_from = c('distance_to_parent','cv_parent','tpm_parent'),names_from = which_parent)
PWdistcheck <- F6_PWwhole %>%
  pivot_wider(values_from = c('distance_to_parent','cv_parent','tpm_parent'),names_from = which_parent)

#Add categories
PCdistcheck$category <- 'Uncategorised'
PCdistcheck$category[PCdistcheck$`distance_to_parent_To C`<0.1 & PCdistcheck$`distance_to_parent_To P`>0.2] <- 'DFO_a'
PCdistcheck$category[PCdistcheck$`distance_to_parent_To C`>0.2 & PCdistcheck$`distance_to_parent_To P`<0.1] <- 'DFO_b'
PCdistcheck$category[PCdistcheck$`distance_to_parent_To C`>0.2 & PCdistcheck$`distance_to_parent_To P`>0.2] <- 'DFB'
PCdistcheck$category[PCdistcheck$`distance_to_parent_To C`<0.2 & PCdistcheck$`distance_to_parent_To P`<0.2] <- 'Conserved'

PWdistcheck$category <- 'Uncategorised'
PWdistcheck$category[PWdistcheck$`distance_to_parent_To W`<0.1 & PWdistcheck$`distance_to_parent_To P`>0.2] <- 'DFO_a'
PWdistcheck$category[PWdistcheck$`distance_to_parent_To W`>0.2 & PWdistcheck$`distance_to_parent_To P`<0.1] <- 'DFO_b'
PWdistcheck$category[PWdistcheck$`distance_to_parent_To W`>0.2 & PWdistcheck$`distance_to_parent_To P`>0.2] <- 'DFB'
PWdistcheck$category[PWdistcheck$`distance_to_parent_To W`<0.2 & PWdistcheck$`distance_to_parent_To P`<0.2] <- 'Conserved'

exPCdistcheck <- PCdistcheck[,c('group_id','genotype','id','category')]
PCsupp_table <- updated_F6_PCwhole %>%
  inner_join(exPCdistcheck)

exPWdistcheck <- PWdistcheck[,c('group_id','genotype','id','category')]
PWsupp_table <- updated_F6_PWwhole %>%
  inner_join(exPWdistcheck)
colnames(PCsupp_table)[1] <- 'f5line'
colnames(PWsupp_table)[1] <- 'f5line'
colnames(PCsupp_table)[18] <- 'bias_distance_to_P'
colnames(PWsupp_table)[18] <- 'bias_distance_to_P'
colnames(PCsupp_table)[19] <- 'bias_distance_to_C'
colnames(PWsupp_table)[19] <- 'bias_distance_to_W'

PCsupp_table <- PCsupp_table[,c('f6line','group_id','A_tpm','B_tpm','D_tpm','A','B','D','cv','associated_homoeologue',
                                'SNP_for_subgenome_collapsed','inherited_genotype_collapsed','bias_distance_to_P','bias_distance_to_C',
                                'category')]
PWsupp_table <- PWsupp_table[,c('f6line','group_id','A_tpm','B_tpm','D_tpm','A','B','D','cv','associated_homoeologue',
                                'SNP_for_subgenome_collapsed','inherited_genotype_collapsed','bias_distance_to_P','bias_distance_to_W',
                                'category')]

PCsupp_table$cv <- PCsupp_table$cv*100
PWsupp_table$cv <- PWsupp_table$cv*100

homologies <- read.csv(file="homoeologs_1_1_1_synt_and_non_synt.csv")
homologies <- homologies[,1:4]
colnames(homologies) <- c('triad_id','A_homoeolog','B_homoeolog','D_homoeolog')

#adjust column names and add gene names
colnames(PCsupp_table) <- c('f5line','triad_id','A_tpm','B_tpm','D_tpm','A_relative_expression','B_relative_expression',
                            'D_relative_expression','cv','homoeolog_associated_expression_with_genotype','SNPs_identified_for_homoeologs_collapsed',
                            'inherited_homoeolog_genotype_from_parent_collapsed','bias_distance_to_Paragon','bias_distance_to_Charger',
                            'bias_distance_category')
colnames(PWsupp_table) <- c('f5line','triad_id','A_tpm','B_tpm','D_tpm','A_relative_expression','B_relative_expression',
                            'D_relative_expression','cv','homoeolog_associated_expression_with_genotype','SNPs_identified_for_homoeologs_collapsed',
                            'inherited_homoeolog_genotype_from_parent_collapsed','bias_distance_to_Paragon','bias_distance_to_Watkins',
                            'bias_distance_category')

PCsupp_table_homoeologs <- PCsupp_table %>%
  inner_join(homologies)
PCsupp_table_homoeologs$inherited_homoeolog_genotype_from_parent_collapsed <- gsub('P1','Paragon',PCsupp_table_homoeologs$inherited_homoeolog_genotype_from_parent_collapsed)
PCsupp_table_homoeologs$inherited_homoeolog_genotype_from_parent_collapsed <- gsub('P2','Charger',PCsupp_table_homoeologs$inherited_homoeolog_genotype_from_parent_collapsed)

PWsupp_table_homoeologs <- PWsupp_table %>%
  inner_join(homologies)
PWsupp_table_homoeologs$inherited_homoeolog_genotype_from_parent_collapsed <- gsub('P1','Paragon',PWsupp_table_homoeologs$inherited_homoeolog_genotype_from_parent_collapsed)
PWsupp_table_homoeologs$inherited_homoeolog_genotype_from_parent_collapsed <- gsub('P2','Watkins',PWsupp_table_homoeologs$inherited_homoeolog_genotype_from_parent_collapsed)

# Now just reorder the columns neatly
PCsupp_table_homoeologs <- PCsupp_table_homoeologs[,c('f5line','triad_id','A_homoeolog','B_homoeolog','D_homoeolog',
                                                      'A_tpm','B_tpm','D_tpm','A_relative_expression','B_relative_expression',
                                                      'D_relative_expression','cv','homoeolog_associated_expression_with_genotype','SNPs_identified_for_homoeologs_collapsed',
                                                      'inherited_homoeolog_genotype_from_parent_collapsed','bias_distance_to_Paragon','bias_distance_to_Charger',
                                                      'bias_distance_category')]
PWsupp_table_homoeologs <- PWsupp_table_homoeologs[,c('f5line','triad_id','A_homoeolog','B_homoeolog','D_homoeolog',
                                                      'A_tpm','B_tpm','D_tpm','A_relative_expression','B_relative_expression',
                                                      'D_relative_expression','cv','homoeolog_associated_expression_with_genotype','SNPs_identified_for_homoeologs_collapsed',
                                                      'inherited_homoeolog_genotype_from_parent_collapsed','bias_distance_to_Paragon','bias_distance_to_Watkins',
                                                      'bias_distance_category')]

write.table(PCsupp_table_homoeologs,'/Users/glombik/work/Documents/my_articles/inheritance_2023/Supplementary_table_2.tsv',quote = F,row.names = F,sep = '\t')
write.table(PWsupp_table_homoeologs,'/Users/glombik/work/Documents/my_articles/inheritance_2023/Supplementary_table_3.tsv',quote = F,row.names = F,sep = '\t')



updated_F6_PWwhole <- read.table('/Users/glombik/work/vcf_retry/updated_F6_PWwhole',header=T,sep = '\t')
updated_F6_PCwhole <- read.table('/Users/glombik/work/vcf_retry/updated_F6_PCwhole',header=T,sep = '\t')
#Now extract only those triads where parental distance is > 0.2

higherdistPC <- read.table('/Users/glombik/work/obj1_reanalysis/highdistPCpar.tsv',header = T)
higherdistPW <- read.table('/Users/glombik/work/obj1_reanalysis/highdistPWpar.tsv',header = T)

higherdist_updated_F6_PCwhole <- updated_F6_PCwhole[updated_F6_PCwhole$group_id %in% higherdistPC$group_id,]
higherdist_updated_F6_PWwhole <- updated_F6_PWwhole[updated_F6_PWwhole$group_id %in% higherdistPW$group_id,]

# Now do plots to manually check...load in parents as well

# Add parental data based on group id in sigtriads
parPC <- read.table("/Users/glombik/work/obj1_reanalysis/HEB_par_PC",header = T)
checkparPC <- parPC[,c("group_id","A_tpm","B_tpm","D_tpm","genotype")]
colnames(checkparPC) <- c("group_id","A_tpm","B_tpm","D_tpm","genotype")
checkparPC$inherited_genotype_collapsed[checkparPC$genotype == "PC_C"] <- "Charger = P2"
checkparPC$inherited_genotype_collapsed[checkparPC$genotype == "PC_P"] <- "Paragon = P1"
checkparPC$associated_homoeologue <- "A,B,D"
checkparPC$SNP_for_subgenome_collapsed <- "A,B,D"
# create same df from F6
checkF6PC <- higherdist_updated_F6_PCwhole[,c("group_id","A_tpm","B_tpm","D_tpm","inherited_genotype_collapsed",
                                              'associated_homoeologue','SNP_for_subgenome_collapsed','genotype')]
checkF6PC$group_id <- as.integer(checkF6PC$group_id)

# Now merge it in

mergedparF6 <- checkparPC[checkparPC$group_id %in% checkF6PC$group_id,] %>%
  rbind(checkF6PC)

# Create a palette that will have spec colour for each genotype (to make it unified)
# First extract all unique GT_filter profiles
GTprofs <- as.data.frame(unique(mergedparF6[,c("inherited_genotype_collapsed")]))
colnames(GTprofs) <- c("inherited_genotype_collapsed")
GTprofs <- GTprofs %>%
  add_row(inherited_genotype_collapsed = 'Watkins = P2')

# palette based on number of profiles
ncolor <- 41
triadpalette <- distinctColorPalette(ncolor)
savepalette <- triadpalette

GTprofs$color <- savepalette[1:length(GTprofs$inherited_genotype_collapsed)]

#merge it in
mergedparF6 <- mergedparF6 %>%
  inner_join(GTprofs)

parentgroup = c("Charger = P2","Paragon = P1")

for (i in unique(mergedparF6$group_id)) {
  v <- ""
  v <- mergedparF6[mergedparF6$group_id==i,]
  d <- ""
  d <- ggtern(v,aes(A_tpm,D_tpm,B_tpm),group=inherited_genotype_collapsed) +
    geom_Tline(Tintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
    geom_Rline(Rintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
    geom_Lline(Lintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
    geom_mask() +
    geom_smooth_tern(data = v %>% filter(!inherited_genotype_collapsed %in% parentgroup),method=lm,fullrange=TRUE,colour='red') +
    geom_point(data = v %>% filter(!inherited_genotype_collapsed %in% parentgroup),aes(fill=inherited_genotype_collapsed),pch=21,colour="black",size=4,alpha=0.8) +
    geom_point(data = v %>% filter(inherited_genotype_collapsed %in% parentgroup),aes(fill=inherited_genotype_collapsed),pch=21,colour="black",size=8,alpha=0.8,show.legend = F) +
    geom_point(data = v %>% filter(!inherited_genotype_collapsed %in% parentgroup),aes(shape=associated_homoeologue),alpha=0) +
    scale_fill_manual(values=unique(v$color[order(v$inherited_genotype_collapsed)])) +
    labs(x="A",y="D",z="B") +
    theme_gray() +
    theme(axis.title = element_text(size = 15),axis.text = element_text(size = 10))
  ggsave(filename = paste0("/Users/glombik/work/vcf_retry/gatk/ggtern/single_homoeolog/PCpassed/group", i,".png"),device = "png",dpi = 400,
         plot = d,width = 8,height = 8,units = "in")
}


# now PW
parPW <- read.table("/Users/glombik/work/obj1_reanalysis/HEB_par_PW",header = T)
checkparPW <- parPW[,c("group_id","A_tpm","B_tpm","D_tpm","genotype")]
colnames(checkparPW) <- c("group_id","A_tpm","B_tpm","D_tpm","genotype")
checkparPW$inherited_genotype_collapsed[checkparPW$genotype == "PW_W"] <- "Watkins = P2"
checkparPW$inherited_genotype_collapsed[checkparPW$genotype == "PW_P"] <- "Paragon = P1"
checkparPW$associated_homoeologue <- "A,B,D"
checkparPW$SNP_for_subgenome_collapsed <- "A,B,D"
# create same df from F6
checkF6PW <- higherdist_updated_F6_PWwhole[,c("group_id","A_tpm","B_tpm","D_tpm","inherited_genotype_collapsed",
                                              'associated_homoeologue','SNP_for_subgenome_collapsed','genotype')]
checkF6PW$group_id <- as.integer(checkF6PW$group_id)

# Now merge it in

mergedparF6PW <- checkparPW[checkparPW$group_id %in% checkF6PW$group_id,] %>%
  rbind(checkF6PW)

# Create a palette that will have spec colour for each genotype (to make it unified)
# First extract all unique GT_filter profiles
GTprofs <- as.data.frame(unique(mergedparF6PW[,c("inherited_genotype_collapsed")]))
colnames(GTprofs) <- c("inherited_genotype_collapsed")
# palette based on number of profiles
# ncolor <- 40
# triadpalette <- distinctColorPalette(ncolor)
# savepalette <- triadpalette

# GTprofs$color <- savepalette[1:length(GTprofs$inherited_genotype_collapsed)]
#merge it in
mergedparF6PW <- mergedparF6PW %>%
  inner_join(GTprofs)

parentgroupPW = c("Watkins = P2","Paragon = P1")

for (i in unique(mergedparF6PW$group_id)) {
  v <- ""
  v <- mergedparF6PW[mergedparF6PW$group_id==i,]
  d <- ""
  d <- ggtern(v,aes(A_tpm,D_tpm,B_tpm),group=inherited_genotype_collapsed) +
    geom_Tline(Tintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
    geom_Rline(Rintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
    geom_Lline(Lintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
    geom_mask() +
    geom_smooth_tern(data = v %>% filter(!inherited_genotype_collapsed %in% parentgroupPW),method=lm,fullrange=TRUE,colour='red') +
    geom_point(data = v %>% filter(!inherited_genotype_collapsed %in% parentgroupPW),aes(fill=inherited_genotype_collapsed),pch=21,colour="black",size=4,alpha=0.8) +
    geom_point(data = v %>% filter(inherited_genotype_collapsed %in% parentgroupPW),aes(fill=inherited_genotype_collapsed),pch=21,colour="black",size=8,alpha=0.8,show.legend = F) +
    geom_point(data = v %>% filter(!inherited_genotype_collapsed %in% parentgroupPW),aes(shape=associated_homoeologue),alpha=0) +
    scale_fill_manual(values=unique(v$color[order(v$inherited_genotype_collapsed)])) +
    labs(x="A",y="D",z="B") +
    theme_gray() +
    theme(axis.title = element_text(size = 15),axis.text = element_text(size = 10))
  ggsave(filename = paste0("/Users/glombik/work/vcf_retry/gatk/ggtern/single_homoeolog/PWpassed/group", i,".png"),device = "png",dpi = 400,
         plot = d,width = 8,height = 8,units = "in")
}


# Check how many triads have how many spec subgenomes associated
numbersPC <- unique(higherdist_updated_F6_PCwhole[,c('group_id','associated_homoeologue')])
table(numbersPC$associated_homoeologue)
# A   A,B A,B,D   A,D     B   B,D     D 
# 115    17     1     9   116     6    44 

numbersPW <- unique(higherdist_updated_F6_PWwhole[,c('group_id','associated_homoeologue')])
table(numbersPW$associated_homoeologue)
# A A,B A,D   B B,D   D 
# 57   4   2  36   2  18


# Lastly, for figure 3, create bias dist plots with a specific triad + triangle plot
distPCexample <- PCdistcheck[PCdistcheck$group_id=='11',]
distPCexample <- distPCexample[,1:2]
distPCexample <- distPCexample %>%
  inner_join(higherdist_updated_F6_PCwhole)

ggplot(distPCexample,aes(`To C`,`To P`)) + 
  geom_point(data=distPCexample,size=5,aes(group=inherited_genotype_collapsed,fill=inherited_genotype_collapsed),pch=21,alpha=0.8,color='black',show.legend = F) +
  scale_color_manual(values = c('#999999','#56B4E9','#009E73','#F0E442','#E69F00')) +
  scale_fill_manual(values = c('#DC4462','#7F9596')) +
  annotate('rect',xmin=0,xmax=0.1,ymin=0.2,ymax=1.2,fill='#56B4E9',alpha=0.15) +
  annotate('rect',xmin=0.2,xmax=1.2,ymin=0,ymax=0.1,fill='#009E73',alpha=0.15) +
  annotate('rect',xmin=0.2,xmax=1.2,ymin=0.2,ymax=1.2,fill='#F0E442',alpha=0.15) +
  annotate('rect',xmin=0,xmax=0.2,ymin=0,ymax=0.2,fill='#E69F00',alpha=0.15) +
  annotate('text',x=0.05,y=0.8,label='DFO_a',size=9,angle=90) +
  annotate('text',x=0.8,y=0.05,label='DFO_b',size=9) +
  annotate('text',x=0.15,y=0.6,label='Uncategorised',size=9,angle=90) +
  annotate('text',x=0.6,y=0.15,label='Uncategorised',size=9) +
  annotate('text',x=0.8,y=0.8,label='DFB',size=9) +
  annotate('text',x=0.1,y=0.1,label='Con',size=9) +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) +
  xlab('Bias distance to Charger parent') +
  ylab('Bias distance to Paragon parent') +
  theme_classic() +
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size=24),
        axis.text = element_text(size = 25))

ggsave('/Users/glombik/work/vcf_retry/gatk/ggtern/single_homoeolog/triad11biasdist.svg',device = 'svg',height = 8,
       width = 8,dpi = 400,units = 'in',plot = last_plot())


##### Try the stats just on CPM data as we thought before!!!
# We expect the opposite = Variation only in the homoeolog that has associated expression with genotype.
# Still stats test p< 0.05 but stdev and variance should be higher for the respective homoeolog!!
# So I need to create a table with triad CPM data combined with inherited genotype data from updatedF6PCwhole data frame +
# significant homoeologs

PCcpm_long <- cpmF6filteredHCcountsPC %>%
  pivot_longer(!gene,names_to = 'f6line',values_to = 'cpm_expression')
colnames(PCcpm_long)[1] <- 'gene.id'
PCcpm_long <- PCcpm_long %>%
  inner_join(lnghomologies,by='gene.id') %>%
  select(c('group_id','genome','cpm_expression','f6line'))

PCcpm_wide <- PCcpm_long %>%
  pivot_wider(names_from = genome,values_from = cpm_expression)
PCcpm_wide$genotype <- substr(PCcpm_wide$f6line,1,5)
PCcpm_wide$genotype <- gsub('_$','',PCcpm_wide$genotype)
colnames(PCcpm_wide)[3:5] <- c('Acpm','Bcpm','Dcpm')

#join with the significant triads
PCcpm_wide <- PCcpm_wide %>%
  inner_join(updated_F6_PCwhole,by=c('group_id','genotype'))
#subset for those, where the parental distance is > 0.2 

PCcpm_wide <- PCcpm_wide[PCcpm_wide$group_id %in% higherdist_updated_F6_PCwhole$group_id,]

PCcpm_angle <- PCcpm_wide %>%
  group_by(group_id) %>%
  summarise(stdev_A=sd(Acpm),stdev_B=sd(Bcpm),stdev_D=sd(Dcpm),variance_A=var(Acpm),variance_B=var(Bcpm),variance_D=var(Dcpm))

PCcpm_angle[is.na(PCcpm_angle)] <- 0

PCcpm_angle$highest_stdev <- apply(PCcpm_angle[,2:4],1,which.max)
PCcpm_angle$highest_stdev <- gsub('1','A',PCcpm_angle$highest_stdev)
PCcpm_angle$highest_stdev <- gsub('2','B',PCcpm_angle$highest_stdev)
PCcpm_angle$highest_stdev <- gsub('3','D',PCcpm_angle$highest_stdev)

PCcpm_angle$highest_var <- apply(PCcpm_angle[,5:7],1,which.max)
PCcpm_angle$highest_var <- gsub('1','A',PCcpm_angle$highest_var)
PCcpm_angle$highest_var <- gsub('2','B',PCcpm_angle$highest_var)
PCcpm_angle$highest_var <- gsub('3','D',PCcpm_angle$highest_var)

PCcpm_angle <- PCcpm_angle %>%
  left_join(select(PCcpm_wide,group_id,associated_homoeologue),by = 'group_id') %>%
  distinct()

#In how many triads is the homoeolog with the highest variation included in the "associated_homoeologue" column

PCcpm_angle$match <- mapply(grepl,PCcpm_angle$highest_var,PCcpm_angle$associated_homoeologue)
# length(unlist(PCcpm_angle[PCcpm_angle$match == TRUE,c('group_id')]))
table(PCcpm_angle$match)
# FALSE  TRUE 
# 124   184 

### Now for PW


PWcpm_long <- cpmF6filteredHCcountsPW %>%
  pivot_longer(!gene,names_to = 'f6line',values_to = 'cpm_expression')
colnames(PWcpm_long)[1] <- 'gene.id'
PWcpm_long <- PWcpm_long %>%
  inner_join(lnghomologies,by='gene.id') %>%
  select(c('group_id','genome','cpm_expression','f6line'))

PWcpm_wide <- PWcpm_long %>%
  pivot_wider(names_from = genome,values_from = cpm_expression)
PWcpm_wide$genotype <- substr(PWcpm_wide$f6line,1,5)
PWcpm_wide$genotype <- gsub('_$','',PWcpm_wide$genotype)
colnames(PWcpm_wide)[3:5] <- c('Acpm','Bcpm','Dcpm')

#join with the significant triads
PWcpm_wide <- PWcpm_wide %>%
  inner_join(updated_F6_PWwhole,by=c('group_id','genotype'))
#subset for those, where the parental distance is > 0.2 

PWcpm_wide <- PWcpm_wide[PWcpm_wide$group_id %in% higherdist_updated_F6_PWwhole$group_id,]

PWcpm_angle <- PWcpm_wide %>%
  group_by(group_id) %>%
  summarise(stdev_A=sd(Acpm),stdev_B=sd(Bcpm),stdev_D=sd(Dcpm),variance_A=var(Acpm),variance_B=var(Bcpm),variance_D=var(Dcpm))

PWcpm_angle[is.na(PWcpm_angle)] <- 0

PWcpm_angle$highest_stdev <- apply(PWcpm_angle[,2:4],1,which.max)
PWcpm_angle$highest_stdev <- gsub('1','A',PWcpm_angle$highest_stdev)
PWcpm_angle$highest_stdev <- gsub('2','B',PWcpm_angle$highest_stdev)
PWcpm_angle$highest_stdev <- gsub('3','D',PWcpm_angle$highest_stdev)

PWcpm_angle$highest_var <- apply(PWcpm_angle[,5:7],1,which.max)
PWcpm_angle$highest_var <- gsub('1','A',PWcpm_angle$highest_var)
PWcpm_angle$highest_var <- gsub('2','B',PWcpm_angle$highest_var)
PWcpm_angle$highest_var <- gsub('3','D',PWcpm_angle$highest_var)

PWcpm_angle <- PWcpm_angle %>%
  left_join(select(PWcpm_wide,group_id,associated_homoeologue),by = 'group_id') %>%
  distinct()

#In how many triads is the homoeolog with the highest variation included in the "associated_homoeologue" column
PWcpm_angle$match <- mapply(grepl,PWcpm_angle$highest_var,PWcpm_angle$associated_homoeologue)
# length(unlist(PWcpm_angle[PWcpm_angle$match == TRUE,c('group_id')]))
table(PWcpm_angle$match)
# FALSE  TRUE 
# 45   74 

# Check the failed ones in both populations. Why do they not follow the pattern?
PCcpm_angle_fail <- PCcpm_angle[PCcpm_angle$match==F,]
PWcpm_angle_fail <- PWcpm_angle[PWcpm_angle$match==F,]

PCcpm_angle_fail <- PCcpm_angle_fail %>%
  left_join(select(PCcpm_wide,group_id,associated_homoeologue,SNP_for_subgenome_collapsed),by = 'group_id') %>%
  distinct()

PWcpm_angle_fail <- PWcpm_angle_fail %>%
  left_join(select(PWcpm_wide,group_id,associated_homoeologue,SNP_for_subgenome_collapsed),by = 'group_id') %>%
  distinct()

PCcpm_angle_fail$snp_geno_no_match <- mapply(grepl,PCcpm_angle_fail$highest_var,PCcpm_angle_fail$SNP_for_subgenome_collapsed)
PCcpm_angle_failT <- PCcpm_angle_fail %>%
  select(group_id,snp_geno_no_match) %>%
  filter(snp_geno_no_match==T) %>%
  distinct()

PWcpm_angle_fail$snp_geno_no_match <- mapply(grepl,PWcpm_angle_fail$highest_var,PWcpm_angle_fail$SNP_for_subgenome_collapsed)
PWcpm_angle_failT <- PWcpm_angle_fail %>%
  select(group_id,snp_geno_no_match) %>%
  filter(snp_geno_no_match==T) %>%
  distinct()


PCcpm_angle_failF <- PCcpm_angle_fail %>%
  select(group_id,snp_geno_no_match,SNP_for_subgenome_collapsed,highest_var) %>%
  filter(snp_geno_no_match==F) %>%
  distinct() %>%
  anti_join(PCcpm_angle_failT,by='group_id')
PWcpm_angle_failF <- PWcpm_angle_fail %>%
  select(group_id,snp_geno_no_match,SNP_for_subgenome_collapsed,highest_var) %>%
  filter(snp_geno_no_match==F) %>%
  distinct() %>%
  anti_join(PWcpm_angle_failT,by='group_id')

uniPCfailF <- PCcpm_angle_failF %>%
  select(group_id,highest_var) %>%
  distinct()
table(uniPCfailF$highest_var)
uniPWfailF <- PWcpm_angle_failF %>%
  select(group_id,highest_var) %>%
  distinct()
table(uniPWfailF$highest_var)

# For figure 3, plot a triad that nicely fits the model and the variation in homoeolog and a second triad where it passes the model
# but fails the variation check showing that we are missing information for other homoeologue/s


goodtriad <- mergedparF6[mergedparF6$group_id=='11',]

ggtern(goodtriad,aes(A_tpm,D_tpm,B_tpm),group=inherited_genotype_collapsed) +
  geom_Tline(Tintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_Rline(Rintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_Lline(Lintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_mask() +
  geom_smooth_tern(data = goodtriad %>% filter(!inherited_genotype_collapsed %in% parentgroup),method=lm,fullrange=TRUE,colour='red') +
  geom_point(data = goodtriad %>% filter(!inherited_genotype_collapsed %in% parentgroup),aes(fill=inherited_genotype_collapsed),pch=21,colour="black",size=4,alpha=0.8) +
  geom_point(data = goodtriad %>% filter(inherited_genotype_collapsed %in% parentgroup),aes(fill=inherited_genotype_collapsed),pch=21,colour="black",size=8,alpha=0.8) +
  geom_point(data = goodtriad %>% filter(!inherited_genotype_collapsed %in% parentgroup),aes(shape=associated_homoeologue),alpha=0) +
  scale_fill_manual(values = c('#fff582ff','#DC4462','#7F9596','#efa9a9ff')) +
  labs(x="A",y="D",z="B") +
  theme_gray() +
  theme(axis.title = element_text(size = 15),axis.text = element_text(size = 10))

ggsave(filename = paste0("/Users/glombik/work/vcf_retry/gatk/ggtern/single_homoeolog/goodtriad11ggtern.svg"),device = "svg",dpi = 400,
       plot = last_plot(),width = 8,height = 8,units = "in")


# plot stdev in associated homoeologues and compare it to the not associated ones in the same triad
head(PCcpm_angle)
PCstdevplot <- PCcpm_angle[PCcpm_angle$highest_stdev==PCcpm_angle$associated_homoeologue,]
PCstdevplots <- PCstdevplot %>%
  pivot_longer(stdev_A:stdev_D,values_to = 'stdev',names_to = 'subgenome_stdev')
PCstdevplots$subgenome_stdev <- gsub('stdev_','',PCstdevplots$subgenome_stdev)
PCstdevplots$logstdev <- log(PCstdevplots$stdev)
PCstdevplots <- PCstdevplots %>%
  mutate(classification = case_when(
    associated_homoeologue == highest_stdev & associated_homoeologue == subgenome_stdev ~ 'associated',
    associated_homoeologue == highest_stdev & associated_homoeologue != subgenome_stdev ~ 'not_associated',
    TRUE ~ NA_character_ # Default case if none of the conditions are met
  ))
PCstdevplots$cross <- 'PxC'
ggplot(PCstdevplots,aes(classification,logstdev,fill=classification)) + geom_boxplot()



head(PWcpm_angle)
PWstdevplot <- PWcpm_angle[PWcpm_angle$highest_stdev==PWcpm_angle$associated_homoeologue,]
PWstdevplots <- PWstdevplot %>%
  pivot_longer(stdev_A:stdev_D,values_to = 'stdev',names_to = 'subgenome_stdev')
PWstdevplots$subgenome_stdev <- gsub('stdev_','',PWstdevplots$subgenome_stdev)
PWstdevplots$logstdev <- log(PWstdevplots$stdev)
PWstdevplots <- PWstdevplots %>%
  mutate(classification = case_when(
    associated_homoeologue == highest_stdev & associated_homoeologue == subgenome_stdev ~ 'associated',
    associated_homoeologue == highest_stdev & associated_homoeologue != subgenome_stdev ~ 'not_associated',
    TRUE ~ NA_character_ # Default case if none of the conditions are met
  ))
PWstdevplots$cross <- 'PxW'

ggplot(PWstdevplots,aes(classification,logstdev,fill=classification)) + geom_boxplot()

bothstdevplot <- rbind(PCstdevplots,PWstdevplots)

#test for significant difference in log(SD) between associated and not associated homoeologs
wilcox_bothstdev <- bothstdevplot[bothstdevplot$logstdev!='-Inf',]
wilcox_bothstdev <- wilcox_bothstdev %>%
  group_by(cross) %>%
  wilcox_test(logstdev ~ classification) %>%
  add_significance() %>%
  add_xy_position(x='cross')
wilcox_bothstdev

ggplot(bothstdevplot,aes(cross,logstdev)) +
  geom_boxplot(aes(color=cross,fill=classification),lwd=2) +
  scale_color_manual(values = c('black','#009128ff')) +
  scale_fill_manual(values = c('white','grey'),labels=c('Associated','Not associated')) +
  scale_x_discrete(labels=c('PxC','PxW')) +
  stat_pvalue_manual(wilcox_bothstdev,label = 'p.signif',tip.length = 0,bracket.nudge.y = 1,label.size = 10) +
  theme_bw() +
  ylab('log(SD) of normalised\nhomoeolog expression (CPM)') +
  xlab('Homoeolog association between inherited genotype\nand expression') +
  theme(axis.title = element_text(size = 32),
        axis.text = element_text(size=26),
        legend.text = element_text(size = 30),
        legend.title = element_blank())
ggsave(plot = last_plot(),filename = 'homoeolog_stdev_stats.svg',device = 'svg',dpi = 400,width = 10,height = 10,units = "in")
