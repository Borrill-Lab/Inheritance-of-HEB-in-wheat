# The goal of this script is to run the eQTL with phenotyping data as a covariate
library(tidyverse)
library(VariantAnnotation)
library(snpsettest)
library(meltr)
library(MatrixEQTL)
library(palettetown)
library(ggplot2)
library(gridExtra)
library(gplots)
library(GGally)
library(network)
library(dplyr)
library(UpSetR)
library(ComplexUpset)
library(lme4)

# First load in PC


# This script will produce SNPs of high quality for eQTL analysis
# Load in data for genotype and depth per sample
all_vcf <- readVcf("all3-final.vcf",param=ScanVcfParam(geno = "GT"))
all_vcfsnv <- isSNV(all_vcf,singleAltOnly=T)
all_vcfGT <- as.data.frame(head(geno(all_vcf))[["GT"]])
all_vcfsnvGT <- subset(all_vcfGT,all_vcfsnv)
all_vcfdp <- readVcf("all3-final.vcf",param=ScanVcfParam(geno = "DP"))
all_vcfdp2 <- as.data.frame(head(geno(all_vcfdp))[["DP"]])
all_vcfsnvdp <- subset(all_vcfdp2,all_vcfsnv)
colnames(all_vcfsnvGT) <- gsub("-","_",colnames(all_vcfsnvGT))
colnames(all_vcfsnvdp) <- gsub("-","_",colnames(all_vcfsnvdp))
colnames(all_vcfsnvdp) <- paste(colnames(all_vcfsnvdp),"DP",sep="_")
all_vcfsnvdp[is.na(all_vcfsnvdp)] <- 0

pcvcfcols <- colnames(all_vcfsnvGT)
# Put them together
all_vcf <- cbind(all_vcfsnvGT,all_vcfsnvdp)
rm(all_vcfdp)
rm(all_vcfGT)
rm(all_vcfdp2)
all_vcf$split <- row.names(all_vcf)
all_vcf <- separate(data = all_vcf,col = split, into = c("chr","REST"),sep = ":")
all_vcf <- separate(data = all_vcf,col = REST, into = c("pos","SNP"),sep = "_")
all_vcf <- separate(data = all_vcf,col = SNP, into = c("REF","ALT"),sep = "/")
rownames(all_vcf) = NULL
# all_vcf$chr <- gsub("chr","",as.character(all_vcf$chr))
all_vcf$chr = substr(all_vcf$chr,4,5)
all_vcf <- all_vcf[all_vcf$chr != "Un",]
# Remove indels by not allowing more then 1 character in the SNP column
all_vcf <- all_vcf %>%
  filter(nchar(REF)<2)
all_vcf$id <- 1:nrow(all_vcf)
all_vcf$pos <- as.numeric(all_vcf$pos)

# # load in the SNP parental-differential data
# diffmap_PC <- read.table("/Users/glombik/work/vcf_retry/diffmap_PC.tab",header = T)
# ex_diffPC <- diffmap_PC %>%
#   dplyr::select(id:ALT)

# DO A NEW FILTERING OF THE VCF FILE based on these criteria
# sites with <3 read depth label as missing data
# remove sites with more than two alleles - lines containing 0/2 1/2 2/2
# remove sites which have more than 50 % GT calls missing
# remove sites which have more than 3 % heterozygote GTs (0/1)
# remove sites with MAF < 0.05
# LASTLY keep only sites which are homozygous in parental samples
rm(all_vcfsnvdp)
rm(all_vcfsnvGT)

depthfilt_all_vcf <- all_vcf
depthfilt_all_vcf[depthfilt_all_vcf=="./."] <- "."
depthfilt_all_vcf[depthfilt_all_vcf==".|."] <- "."
depthfilt_all_vcf[depthfilt_all_vcf=="0|0"] <- "0/0"
depthfilt_all_vcf[depthfilt_all_vcf=="1|1"] <- "1/1"
depthfilt_all_vcf[depthfilt_all_vcf=="0|1"] <- "0/1"

# Labeling sites with read depth <3 as missing data
for (i in pcvcfcols) {
  snPCol <- i
  dPCol <- paste0(i,'_DP')
  depthfilt_all_vcf[[snPCol]][depthfilt_all_vcf[[dPCol]]<3] <- "."
}

write.table(depthfilt_all_vcf,file = 'PCdepthfilt_all_vcf',quote = F,row.names = T,sep = '\t')
depthfilt_all_vcf <- read.table('PCdepthfilt_all_vcf',header = T)

# Removing sites that are not homozygous in parents
parentfilt_all_vcf <- depthfilt_all_vcf
parentfilt_all_vcf <- parentfilt_all_vcf[!(parentfilt_all_vcf$PC_P_3=="0/1" | parentfilt_all_vcf$PC_C_3=="0/1"),]
parentfilt_all_vcf <- parentfilt_all_vcf[!(parentfilt_all_vcf$PC_P_3=="." | parentfilt_all_vcf$PC_C_3=="."),]
parentfilt_all_vcf <- parentfilt_all_vcf[!(parentfilt_all_vcf$PC_P_3==parentfilt_all_vcf$PC_C_3),]
#57 700
# Removing sites with more than two alleles - lines containing 0/2 1/2 2/2
sitefilt_all_vcf <- parentfilt_all_vcf
sitefilt_all_vcf <- sitefilt_all_vcf[rowSums(sitefilt_all_vcf[1:162] == "0/2") == 0, ]
sitefilt_all_vcf <- sitefilt_all_vcf[rowSums(sitefilt_all_vcf[1:162] == "1/2") == 0, ]
sitefilt_all_vcf <- sitefilt_all_vcf[rowSums(sitefilt_all_vcf[1:162] == "2/2") == 0, ]
#57 700
# Removing sites with more than 50 % calls missing
sitemissfilt_all_vcf <- sitefilt_all_vcf
sitemissfilt_all_vcf$sitemiss <- rowSums(sitemissfilt_all_vcf[1:160] == ".")
sitemissfilt_all_vcf <- sitemissfilt_all_vcf[sitemissfilt_all_vcf$sitemiss<((160/100)*50),1:329]
# 46 229 sites left

# Removing sites with more than 3 % calls heterozygous
hetfilt_all_vcf <- sitemissfilt_all_vcf
# rm(depthfilt_all_vcf)
hetfilt_all_vcf$het <- rowSums(hetfilt_all_vcf[1:160]=="0/1")
hetfilt_all_vcf <- hetfilt_all_vcf[hetfilt_all_vcf$het<=((160/100)*3),1:329]
# 11 970 sites left
# # Removing sites with MAF < 0.05
# # Formula: 1-((nAA + 0.5 * nAB)/(nAA + nAB + nBB))
maffilt_all_vcf <- hetfilt_all_vcf
maffilt_all_vcf$AA <- rowSums(maffilt_all_vcf[1:160] == "0/0")
maffilt_all_vcf$AB <- rowSums(maffilt_all_vcf[1:160] == "0/1")
maffilt_all_vcf$BB <- rowSums(maffilt_all_vcf[1:160] == "1/1")
maffilt_all_vcf$MAF <- 1-((maffilt_all_vcf$AA + 0.5 * maffilt_all_vcf$AB)/(maffilt_all_vcf$AA + maffilt_all_vcf$AB + maffilt_all_vcf$BB))
maffilt_all_vcf <- maffilt_all_vcf %>%
  filter(rowSums(dplyr::select(.,c('AA','AB','BB'))>=8)>=2)
maffilt_all_vcf <- maffilt_all_vcf[maffilt_all_vcf$MAF>0.05 & maffilt_all_vcf$MAF<0.95,1:329]
#11 283 sites left

# And connect it to the annotation file to prepare the base map for viewing
TraeCSgff <- read.table("/Users/glombik/work/vcf_retry/geneIWGSC_v1.1_HC_20170706.gff3sorted")
head(TraeCSgff)
colnames(TraeCSgff) <- c("chr","start","end","gene.id","strand")
TraeCSgff$chr = substr(TraeCSgff$chr,4,5)
TraeCSgff <- TraeCSgff[TraeCSgff$chr != "Un",]

# Now map the SNPs to the annotation
map_all <- map_snp_to_gene(maffilt_all_vcf,TraeCSgff,extend_start = 0,extend_end = 0)
map_all <- as.data.frame(map_all[["map"]])


map_all$id <- as.numeric(map_all$id)
map_all <- map_all %>%
  dplyr::select(!(gene.adj.start:gene.adj.end)) %>%
  inner_join(sitemissfilt_all_vcf)
# Remove NAs - keep only assigned SNPs
map_all <- map_all[complete.cases(map_all),]
# How many are gene unique?
uniqids <- map_all %>%
  dplyr::count(id)
uniqids_one <- uniqids[uniqids$n == 1,]

# where does the SNP belong - take only SNPs that do not overlap between annotated genes
uni_map_all <- map_all %>%
  inner_join(uniqids_one)
uni_map_all$id <- as.numeric(uni_map_all$id)
map_all <- uni_map_all %>%
  dplyr::select(!(n)) %>%
  inner_join(sitemissfilt_all_vcf)

eQTL_PC <- map_all
#het filt 3% = 7957 SNPS left

# do PCA to see if they cluster by SNPs
PCa_PC_snp <- eQTL_PC %>%
  dplyr::select(c('id',starts_with('PC'))) %>%
  dplyr::select(!ends_with('DP'))
rownames(PCa_PC_snp) <- PCa_PC_snp$id

ids_PCa <- as.data.frame(rownames(PCa_PC_snp))
colnames(ids_PCa) <- 'id'
ids_PCa$id <- as.integer(ids_PCa$id)


PCa_PC_snp <- PCa_PC_snp %>%
  dplyr::select(!id)

# create a heatmap to test check the missing data
PCa_PC_snp <- PCa_PC_snp %>%
  mutate_all(~ifelse(.== PC_P_3,2,.)) %>%
  mutate_all(~ifelse(.==PC_C_3,4,.))



PCa_PC_snp[PCa_PC_snp=="."] <- 0
# PCa_PW_snp[PCa_PW_snp=="0/0"] <- as.numeric(3)
PCa_PC_snp[PCa_PC_snp=="0/1"] <- as.numeric(3)
# PCa_PW_snp[PCa_PW_snp=="1/1"] <- as.numeric(-1)
PCa_PC_snp <- PCa_PC_snp %>%
  mutate_all(as.numeric)

library(gplots)

filteex <- PCa_PC_snp

ggfilteex <- PCa_PC_snp
ggfilteex$snpid <- rownames(ggfilteex)
ggfilteex <- ggfilteex %>%
  pivot_longer(!snpid,values_to = 'GT',names_to = 'sample')

chrsplits <- c(87786,188546,239040,352312,466185,532619,608940,728570,774098,848184,900835,967782,
               1002473,1101982,1146098,1229013,1330746,1371916,1465673,1549152)
ggfilteex$GT <- as.factor(ggfilteex$GT)
# ggfilteex$snpid <- as.numeric(ggfilteex$snpid)
ggplot(ggfilteex,aes(snpid,sample,fill=GT)) + geom_tile() +
  scale_fill_manual(values = c('black','yellow','orange','red')) +
  # geom_vline(xintercept = as.numeric(factor(chrsplits,levels = ggfilteex$snpid)),color='black',linetype='dashed',linewidth=17) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 1))
ggsave("heatcheckPC_recoded.png",device = "png",dpi = 400,plot = last_plot(),width = 20,height = 20,units = "in")


transposed_data <- t(filteex)

pdf('heatcheckPCrecoded_clust.pdf')
heatmap.2(
  transposed_data,
  col = colorRampPalette(c('black','yellow','orange','red'))(100), # Color palette
  scale = "none",  # Scale rows (genes)
  main = "Gene Expression Heatmap",
  trace = 'none',
  cexRow = 0.3,
  cexCol = 0.3,
  dendrogram = "row",
  Colv=F
)
dev.off()



PCa_PC_snp[PCa_PC_snp=="."] <- NA
PCa_PC_snp[PCa_PC_snp=="0/0"] <- as.numeric(1)
PCa_PC_snp[PCa_PC_snp=="0/1"] <- as.numeric(2)
PCa_PC_snp[PCa_PC_snp=="1/1"] <- as.numeric(3)
PCa_PC_snp <- PCa_PC_snp[complete.cases(PCa_PC_snp),]
PCa_PC_snp <- PCa_PC_snp %>%
  mutate_all(as.numeric)



zsnp <- scale(PCa_PC_snp,center = T,scale = T)
tsnp <- t(zsnp)
snp_pca <- prcomp(tsnp)
snp_out <- as.data.frame(snp_pca$x)
snp_out$group <- sapply(strsplit(as.character(row.names(tsnp)),"_"),"[[",1)
head(snp_out)
percentage <- round(snp_pca$sdev / sum(snp_pca$sdev) * 100, 2)
percentage <- paste( colnames(snp_out), "(", paste( as.character(percentage), "%", ")", sep="") )
parents <- c('PC_C_3','PC_P_3','PW_P_3','PW_W_3')

## Also take 20 PCs as a covariate for eQTL with all samples?
covar_pc1to3 <- as.data.frame(t(snp_out[,1:20]))
colnames(covar_pc1to3) <- rownames(snp_out)
covar_pc1to3$id <- rownames(covar_pc1to3)
covar_pc1to3 <- covar_pc1to3 %>%
  relocate(id,.before = PC_1_3)
covar_pc1to3 <- covar_pc1to3 %>%
  dplyr::select(!matches(c('PC_P','PC_C')))
colnames(covar_pc1to3) <- gsub("_$", "", colnames(covar_pc1to3))
covar_pc1to3_order <- order(colnames(covar_pc1to3))
covar_pc1to3 <- covar_pc1to3[,covar_pc1to3_order]
write.table(covar_pc1to3,file = 'newwayPCC3covar_pc1to20',row.names = F,quote = F,sep = '\t')
##
p <- ggplot(snp_out,aes(x=PC1,y=PC2,color=group)) +
  geom_point(size=4.5,alpha=0.5) +
  geom_text(aes(label = rownames(snp_out))) +
  geom_text(data=snp_out %>% filter(rownames(snp_out) %in% parents),aes(label=rownames(snp_out %>% filter(rownames(snp_out) %in% parents))),colour='black') +
  xlab(percentage[1]) +
  ylab(percentage[2])
p

#remove the parents
eQTL_PC <- map_all %>%
  dplyr::select(!matches(c('PC_P','PC_C')))

write.table(eQTL_PC,file = 'gatknewway_filteQTLPCC3',row.names = F,quote = F,sep = '\t')





################### RUN THE EQTL
filteQTLPC <- read.table('gatknewway_filteQTLPCC3',header = T)
colnames(filteQTLPC) <- gsub("_$", "", colnames(filteQTLPC))

gt_PC <- filteQTLPC %>%
  dplyr::select(c('id',starts_with('PC'))) %>%
  dplyr::select(!ends_with('DP'))
gt_PC[gt_PC=="0/0"] <- as.numeric(0)
gt_PC[gt_PC=="0/1"] <- as.numeric(1)
gt_PC[gt_PC=="1/1"] <- as.numeric(2)
gt_PC[gt_PC=="."] <- NA
gtPC_order <- order(colnames(gt_PC))
gt_PC <- gt_PC[,gtPC_order]
colnames(gt_PC) <- gsub("_$", "", colnames(gt_PC))
write.table(gt_PC,file = 'gatknewway_gt_PCC3',sep = '\t',row.names = F,quote = F)

# Gene position data
gene_loc <- filteQTLPC %>%
  dplyr::select(c('gene.id','chr','gene.start','gene.end'))
write.table(gene_loc,file = 'newway_gene_locPCC3',sep = '\t',row.names = F,quote = F)

#SNP position data
snp_loc <- filteQTLPC %>%
  dplyr::select(c('id','chr','pos'))
write.table(snp_loc,file = 'newway_snp_locPCC3',sep = '\t',row.names = F,quote = F)







##### Now do the same for PW population

# Load in data for genotype and depth per sample
all_vcf <- readVcf("gatk/PW/all3.vcf",param=ScanVcfParam(geno = "GT"))
all_vcfsnv <- isSNV(all_vcf,singleAltOnly=T)
all_vcfGT <- as.data.frame(head(geno(all_vcf))[["GT"]])
all_vcfsnvGT <- subset(all_vcfGT,all_vcfsnv)
all_vcfdp <- readVcf("gatk/PW/all3.vcf",param=ScanVcfParam(geno = "DP"))
all_vcfdp2 <- as.data.frame(head(geno(all_vcfdp))[["DP"]])
all_vcfsnvdp <- subset(all_vcfdp2,all_vcfsnv)
colnames(all_vcfsnvGT) <- gsub("-","_",colnames(all_vcfsnvGT))
colnames(all_vcfsnvdp) <- gsub("-","_",colnames(all_vcfsnvdp))
colnames(all_vcfsnvdp) <- paste(colnames(all_vcfsnvdp),"DP",sep="_")
all_vcfsnvdp[is.na(all_vcfsnvdp)] <- 0

pcvcfcols <- colnames(all_vcfsnvGT)
# Put them together
all_vcf <- cbind(all_vcfsnvGT,all_vcfsnvdp)

all_vcf$split <- row.names(all_vcf)
all_vcf <- separate(data = all_vcf,col = split, into = c("chr","REST"),sep = ":")
all_vcf <- separate(data = all_vcf,col = REST, into = c("pos","SNP"),sep = "_")
all_vcf <- separate(data = all_vcf,col = SNP, into = c("REF","ALT"),sep = "/")
rownames(all_vcf) = NULL
# all_vcf$chr <- gsub("chr","",as.character(all_vcf$chr))
all_vcf$chr = substr(all_vcf$chr,4,5)
all_vcf <- all_vcf[all_vcf$chr != "Un",]
# Remove indels by not allowing more then 1 character in the SNP column
all_vcf <- all_vcf %>%
  filter(nchar(REF)<2)
all_vcf$id <- 1:nrow(all_vcf)
all_vcf$pos <- as.numeric(all_vcf$pos)

# # load in the SNP parental-differential data
# diffmap_PC <- read.table("/Users/glombik/work/vcf_retry/diffmap_PC.tab",header = T)
# ex_diffPC <- diffmap_PC %>%
#   dplyr::select(id:ALT)

# DO A NEW FILTERING OF THE VCF FILE based on these criteria
# sites with <3 read depth label as missing data
# remove sites with more than two alleles - lines containing 0/2 1/2 2/2
# remove sites which have more than 50 % GT calls missing
# remove sites which have more than 3 % heterozygote GTs (0/1)
# remove sites with MAF < 0.05
# LASTLY keep only sites which are homozygous in parental samples
rm(all_vcfdp)
rm(all_vcfdp2)
rm(all_vcfGT)
rm(all_vcfsnvdp)
rm(all_vcfsnvGT)


# Filter for batch1 already here
# But dont forget to extract the first columns with id,chr,...!!!
pwallbatchdf <- read.table('pwallbatchdf',header=F)
batch1ids <- pwallbatchdf[1,pwallbatchdf[2,]==1]
batch1ids <- unlist(batch1ids)
importantcols <- c('chr','pos','REF','ALT','id')
batch1ids <- c(batch1ids,importantcols)
b1all_vcf <- all_vcf %>%
  dplyr::select(matches(batch1ids))
column_names <- colnames(b1all_vcf)
ordered_column_names <- c('chr','pos','REF','ALT','id',column_names[!grepl("DP", column_names)], column_names[grepl("DP", column_names)])
b1all_vcf <- b1all_vcf[, ordered_column_names]
b1all_vcf <- b1all_vcf %>%
  dplyr::select(!c(chr.1,pos.1,REF.1,ALT.1,id.1))
row1allvcf <- b1all_vcf[1,]
rm(all_vcf)


depthfilt_all_vcf <- b1all_vcf

depthfilt_all_vcf[depthfilt_all_vcf=="./."] <- "."
depthfilt_all_vcf[depthfilt_all_vcf==".|."] <- "."
depthfilt_all_vcf[depthfilt_all_vcf=="0|0"] <- "0/0"
depthfilt_all_vcf[depthfilt_all_vcf=="1|1"] <- "1/1"
depthfilt_all_vcf[depthfilt_all_vcf=="0|1"] <- "0/1"


# Labeling sites with read depth <3 as missing data
for (i in column_names) {
  snPCol <- i
  dPCol <- paste0(i,'_DP')
  depthfilt_all_vcf[[snPCol]][depthfilt_all_vcf[[dPCol]]<3] <- "."
}

write.table(depthfilt_all_vcf,file = 'PWdepthfilt_all_vcf_b1',quote = F,row.names = T,sep = '\t')
depthfilt_all_vcf <- read.table('PWdepthfilt_all_vcf_b1',header = T)

# Removing sites that are not homozygous in parents
parentfilt_all_vcf <- depthfilt_all_vcf
parentfilt_all_vcf <- parentfilt_all_vcf[!(parentfilt_all_vcf$PW_P_3_=="0/1" | parentfilt_all_vcf$PW_W_3_=="0/1"),]
parentfilt_all_vcf <- parentfilt_all_vcf[!(parentfilt_all_vcf$PW_P_3_=="." | parentfilt_all_vcf$PW_W_3_=="."),]
parentfilt_all_vcf <- parentfilt_all_vcf[!(parentfilt_all_vcf$PW_P_3_==parentfilt_all_vcf$PW_W_3_),]
#30 187
# Removing sites with more than two alleles - lines containing 0/2 1/2 2/2
sitefilt_all_vcf <- parentfilt_all_vcf
sitefilt_all_vcf <- sitefilt_all_vcf[rowSums(sitefilt_all_vcf[6:57] == "0/2") == 0, ]
sitefilt_all_vcf <- sitefilt_all_vcf[rowSums(sitefilt_all_vcf[6:57] == "1/2") == 0, ]
sitefilt_all_vcf <- sitefilt_all_vcf[rowSums(sitefilt_all_vcf[6:57] == "2/2") == 0, ]
#30 187
# Removing sites with more than 50 % calls missing
sitemissfilt_all_vcf <- sitefilt_all_vcf
sitemissfilt_all_vcf$sitemiss <- rowSums(sitemissfilt_all_vcf[6:55] == ".")
sitemissfilt_all_vcf <- sitemissfilt_all_vcf[sitemissfilt_all_vcf$sitemiss<((50/100)*50),1:109]
# 24 538 sites left

# Removing sites with more than 3 % calls heterozygous
hetfilt_all_vcf <- sitemissfilt_all_vcf
rm(depthfilt_all_vcf)
hetfilt_all_vcf$het <- rowSums(hetfilt_all_vcf[6:55]=="0/1")
hetfilt_all_vcf <- hetfilt_all_vcf[hetfilt_all_vcf$het<=((50/100)*3),1:109]
# 8 040
# # Removing sites with MAF < 0.05
# # Formula: 1-((nAA + 0.5 * nAB)/(nAA + nAB + nBB))
maffilt_all_vcf <- hetfilt_all_vcf
maffilt_all_vcf$AA <- rowSums(maffilt_all_vcf[6:55] == "0/0")
maffilt_all_vcf$AB <- rowSums(maffilt_all_vcf[6:55] == "0/1")
maffilt_all_vcf$BB <- rowSums(maffilt_all_vcf[6:55] == "1/1")
maffilt_all_vcf$MAF <- 1-((maffilt_all_vcf$AA + 0.5 * maffilt_all_vcf$AB)/(maffilt_all_vcf$AA + maffilt_all_vcf$AB + maffilt_all_vcf$BB))
maffilt_all_vcf <- maffilt_all_vcf %>%
  filter(rowSums(dplyr::select(.,c('AA','AB','BB'))>=3)>=2)
maffilt_all_vcf <- maffilt_all_vcf[maffilt_all_vcf$MAF>0.05 & maffilt_all_vcf$MAF<0.95,1:109]
#7 868
maffilt_all_vcf <- maffilt_all_vcf[complete.cases(maffilt_all_vcf),]
maffilt_all_vcf$id <- as.numeric(maffilt_all_vcf$id)
maffilt_all_vcf$pos <- as.numeric(maffilt_all_vcf$pos)

# And connect it to the annotation file to prepare the base map for viewing
TraeCSgff <- read.table("/Users/glombik/work/vcf_retry/geneIWGSC_v1.1_HC_20170706.gff3sorted")
head(TraeCSgff)
colnames(TraeCSgff) <- c("chr","start","end","gene.id","strand")
TraeCSgff$chr = substr(TraeCSgff$chr,4,5)
TraeCSgff <- TraeCSgff[TraeCSgff$chr != "Un",]

# Now map the SNPs to the annotation
map_all <- map_snp_to_gene(maffilt_all_vcf,TraeCSgff,extend_start = 0,extend_end = 0)
map_all <- as.data.frame(map_all[["map"]])


map_all$id <- as.numeric(map_all$id)
map_all <- map_all %>%
  dplyr::select(!(gene.adj.start:gene.adj.end)) %>%
  inner_join(maffilt_all_vcf)
# Remove NAs - keep only assigned SNPs
map_all <- map_all[complete.cases(map_all),]
# How many are gene unique?
uniqids <- map_all %>%
  dplyr::count(id)
uniqids_one <- uniqids[uniqids$n == 1,]

# where does the SNP belong - take only SNPs that do not overlap between annotated genes
uni_map_all <- map_all %>%
  inner_join(uniqids_one)
uni_map_all$id <- as.numeric(uni_map_all$id)
map_all <- uni_map_all %>%
  dplyr::select(!(n)) %>%
  inner_join(maffilt_all_vcf)

#5 538
eQTL_PW <- map_all


# do PCA to see if they cluster by SNPs
PCa_PW_snp <- eQTL_PW %>%
  dplyr::select(c('id',starts_with('PW'))) %>%
  dplyr::select(!ends_with('DP'))
rownames(PCa_PW_snp) <- PCa_PW_snp$id

ids_PCa <- as.data.frame(rownames(PCa_PW_snp))
colnames(ids_PCa) <- 'id'
ids_PCa$id <- as.integer(ids_PCa$id)


PCa_PW_snp <- PCa_PW_snp %>%
  dplyr::select(!id)


# create a heatmap to test check the missing data
PCa_PW_snp <- PCa_PW_snp %>%
  mutate_all(~ifelse(.== PW_P_3_,2,.)) %>%
  mutate_all(~ifelse(.==PW_W_3_,4,.))



PCa_PW_snp[PCa_PW_snp=="."] <- 0
# PCa_PW_snp[PCa_PW_snp=="0/0"] <- as.numeric(3)
PCa_PW_snp[PCa_PW_snp=="0/1"] <- as.numeric(3)
# PCa_PW_snp[PCa_PW_snp=="1/1"] <- as.numeric(-1)
PCa_PW_snp <- PCa_PW_snp %>%
  mutate_all(as.numeric)

library(gplots)

filteex <- PCa_PW_snp

ggfilteex <- PCa_PW_snp
ggfilteex$snpid <- rownames(ggfilteex)
ggfilteex <- ggfilteex %>%
  pivot_longer(!snpid,values_to = 'GT',names_to = 'sample')

chrsplits <- c(87786,188546,239040,352312,466185,532619,608940,728570,774098,848184,900835,967782,
               1002473,1101982,1146098,1229013,1330746,1371916,1465673,1549152)
ggfilteex$GT <- as.factor(ggfilteex$GT)
# ggfilteex$snpid <- as.numeric(ggfilteex$snpid)
ggplot(ggfilteex,aes(snpid,sample,fill=GT)) + geom_tile() +
  scale_fill_manual(values = c('black','yellow','orange','red')) +
  # geom_vline(xintercept = as.numeric(factor(chrsplits,levels = ggfilteex$snpid)),color='black',linetype='dashed',linewidth=17) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 1))
ggsave("heatcheckPW_recoded_b1.png",device = "png",dpi = 400,plot = last_plot(),width = 20,height = 20,units = "in")


transposed_data <- t(filteex)
heatmap.2(
  transposed_data,
  col = colorRampPalette(c('black','yellow','orange','red'))(100), # Color palette
  scale = "none",  # Scale rows (genes)
  main = "Gene Expression Heatmap",
  trace = 'none',
  cexRow = 0.3,
  cexCol = 0.3,
  dendrogram = "row",
  Colv=F
)



# Now the PCA
PCa_PW_snp <- eQTL_PW %>%
  dplyr::select(c('id',starts_with('PW'))) %>%
  dplyr::select(!ends_with('DP'))
rownames(PCa_PW_snp) <- PCa_PW_snp$id

ids_PCa <- as.data.frame(rownames(PCa_PW_snp))
colnames(ids_PCa) <- 'id'
ids_PCa$id <- as.integer(ids_PCa$id)


PCa_PW_snp <- PCa_PW_snp %>%
  dplyr::select(!id)



PCa_PW_snp[PCa_PW_snp=="."] <- NA
PCa_PW_snp[PCa_PW_snp=="0/0"] <- as.numeric(1)
PCa_PW_snp[PCa_PW_snp=="0/1"] <- as.numeric(2)
PCa_PW_snp[PCa_PW_snp=="1/1"] <- as.numeric(3)
PCa_PW_snp <- PCa_PW_snp %>%
  mutate_all(as.numeric)


dPCa_PW_snp <- PCa_PW_snp[complete.cases(PCa_PW_snp),]
dPCa_PW_snp <- dPCa_PW_snp %>%
  mutate_all(as.numeric)


zsnp <- scale(dPCa_PW_snp,center = T,scale = T)
tsnp <- t(zsnp)

snp_pca <- prcomp(tsnp)
snp_out <- as.data.frame(snp_pca$x)
snp_out$group <- sapply(strsplit(as.character(row.names(tsnp)),"_"),"[[",1)
head(snp_out)
percentage <- round(snp_pca$sdev / sum(snp_pca$sdev) * 100, 2)
percentage <- paste( colnames(snp_out), "(", paste( as.character(percentage), "%", ")", sep="") )
parents <- c('PC_C_2','PC_P_3','PW_P_3_','PW_W_3_')

p <- ggplot(snp_out,aes(x=PC1,y=PC2,color=group)) +
  geom_point(size=4.5,alpha=0.5) +
  geom_text(aes(label = rownames(snp_out))) +
  geom_text(data=snp_out %>% filter(rownames(snp_out) %in% parents),aes(label=rownames(snp_out %>% filter(rownames(snp_out) %in% parents))),colour='black') +
  xlab(percentage[1]) +
  ylab(percentage[2])

p



## Also take 20 PCs as a covariate for eQTL with all samples
covar_pw1to3 <- as.data.frame(t(snp_out[,1:20]))
colnames(covar_pw1to3) <- rownames(snp_out)
covar_pw1to3$id <- rownames(covar_pw1to3)
covar_pw1to3 <- covar_pw1to3 %>%
  relocate(id,.before = PW_10_2_)
colnames(covar_pw1to3) <- gsub("__$", "", colnames(covar_pw1to3))
colnames(covar_pw1to3) <- gsub("_$", "", colnames(covar_pw1to3))
covar_pw1to3 <- covar_pw1to3 %>%
  dplyr::select(!matches(c('PW_P','PW_W')))
covar_pw1to3_order <- order(colnames(covar_pw1to3))
covar_pw1to3 <- covar_pw1to3[,covar_pw1to3_order]
write.table(covar_pw1to3,file = 'gatkcovar_pw1to20_batch1',row.names = F,quote = F,sep = '\t')


#remove the parents
eQTL_PW <- map_all %>%
  dplyr::select(!matches(c('PW_P','PW_W')))

write.table(eQTL_PW,file = 'gatknewway_filteQTLPW_batch1',row.names = F,quote = F,sep = '\t')


filteQTLPW <- read.table('gatknewway_filteQTLPW_batch1',header = T)

gt_PW <- filteQTLPW %>%
  dplyr::select(c('id',starts_with('PW'))) %>%
  dplyr::select(!ends_with('DP'))
gt_PW[gt_PW=="0/0"] <- as.numeric(0)
gt_PW[gt_PW=="0/1"] <- as.numeric(1)
gt_PW[gt_PW=="1/1"] <- as.numeric(2)
gt_PW[gt_PW=="."] <- NA
gtPW_order <- order(colnames(gt_PW))
gt_PW <- gt_PW[,gtPW_order]
colnames(gt_PW) <- gsub("__$", "", colnames(gt_PW))
colnames(gt_PW) <- gsub("_$", "", colnames(gt_PW))
# gt_PW <- gt_PW[,which(colnames(gt_PW) %in% covar_samples_left)]
# gt_PW <- gt_PW[,which(colnames(gt_PW) %in% covar_samples_left_minusPar)]
write.table(gt_PW,file = 'gatknewway_gt_PW_batch1',sep = '\t',row.names = F,quote = F)

# Gene position data
gene_loc <- filteQTLPW %>%
  dplyr::select(c('gene.id','chr','gene.start','gene.end'))
write.table(gene_loc,file = 'gatknewway_gene_locPW_batch1',sep = '\t',row.names = F,quote = F)

#SNP position data
snp_loc <- filteQTLPW %>%
  dplyr::select(c('id','chr','pos'))
write.table(snp_loc,file = 'gatknewway_snp_locPW_batch1',sep = '\t',row.names = F,quote = F)




############## run separate model with covariates first and then Matrix eQTL
# ANd change covariate from phenotype values to tray column!


# First the new phenotype file
filteQTLPC <- read.table('Detection_of_eQTL/gatknewway_filteQTLPCC3',header = T)



rev_gatknewwayPC <- read.table('Detection_of_eQTL/gatknewway_gt_PCC3',header = T)

rev_pcallbatchdf <- read.table('Detection_of_eQTL/new_prep_data/f6allbatchdf',header = T) 
rev_pcallbatchdf <- rev_pcallbatchdf[rev_pcallbatchdf$batch==1 & rev_pcallbatchdf$cross==1,]
rev_pcallbatchdf_cols <- unlist(rownames(rev_pcallbatchdf))
rev_pcallbatchdf_cols <- c('id',rev_pcallbatchdf_cols)
rev_gatknewwayPC <- rev_gatknewwayPC[,names(rev_gatknewwayPC) %in% rev_pcallbatchdf_cols]
write.table(rev_gatknewwayPC,'Detection_of_eQTL/new_prep_data/rev_gatknewwayPC_excel',quote = F,row.names = F,sep = '\t')

rev_gatknewwayPC_batch2 <- read.table('Detection_of_eQTL/gatknewway_gt_PCC3',header = T)
rev_pcallbatchdf_cols_batch2 <- read.table('Detection_of_eQTL/new_prep_data/f6allbatchdf',header = T) 
rev_pcallbatchdf_cols_batch2 <- rev_pcallbatchdf_cols_batch2[rev_pcallbatchdf_cols_batch2$batch==2 & rev_pcallbatchdf_cols_batch2$cross==1,]
rev_pcallbatchdf_cols_batch2 <- unlist(rownames(rev_pcallbatchdf_cols_batch2))
rev_pcallbatchdf_cols_batch2 <- c('id',rev_pcallbatchdf_cols_batch2)

rev_gatknewwayPC_batch2 <- rev_gatknewwayPC_batch2[,names(rev_gatknewwayPC_batch2) %in% rev_pcallbatchdf_cols_batch2]
# Load in new GT data and create principal components for covariate file
corr_rev_gatknewwayPC <- read.table('Detection_of_eQTL/new_prep_data/correct_rev_gatknewwayPC_excel',header = T)


#merge with batch 2 data
corr_rev_both_batchesPC <- corr_rev_gatknewwayPC %>%
  inner_join(rev_gatknewwayPC_batch2,by='id')

#generate principal components
PCA_rev_gt <- corr_rev_both_batchesPC

PCA_rev_gt[PCA_rev_gt=="."] <- NA
PCA_rev_gt <- PCA_rev_gt[complete.cases(PCA_rev_gt),]
PCA_rev_gt <- PCA_rev_gt %>%
  mutate_all(as.numeric)

zsnp <- scale(PCA_rev_gt,center = T,scale = T)
tsnp <- t(zsnp)
snp_pca <- prcomp(tsnp)
snp_out <- as.data.frame(snp_pca$x)
snp_out$group <- sapply(strsplit(as.character(row.names(tsnp)),"_"),"[[",1)
head(snp_out)

## Also take 20 PCs as a covariate for eQTL with all samples?
rev_covar_PC20 <- as.data.frame(t(snp_out[,1:20]))
colnames(rev_covar_PC20) <- rownames(snp_out)
rev_covar_PC20$id <- rownames(rev_covar_PC20)
#column order for all files now
revcovar_pc1to3_order <- order(colnames(rev_covar_PC20))
rev_covar_PC20 <- rev_covar_PC20[,revcovar_pc1to3_order]


#order the GT file as well
corr_rev_both_batchesPC <- corr_rev_both_batchesPC[,revcovar_pc1to3_order]
write.table(corr_rev_both_batchesPC,'Detection_of_eQTL/new_prep_data/ready_rev_gt_PC_both_batches',quote = F,row.names = F,sep = '\t')

#now add batch info to the covariate file

batchinfo <- read.table('Detection_of_eQTL/new_prep_data/f6allbatchdf',header = T)
batchinfo <- batchinfo[batchinfo$cross==1,]
batchinfo <- as.data.frame(t(batchinfo))
batchinfo$id <- rownames(batchinfo)
batchinfo <- batchinfo %>%
  select(id,everything())
rownames(batchinfo) <- NULL
batchinfo <- batchinfo %>%
  filter(!row_number() %in% c(2))

select_pcs <- colnames(rev_covar_PC20)

batchinfo <- batchinfo[,match(select_pcs, colnames(batchinfo))]

rev_covar_PC20 <- rbind(rev_covar_PC20,batchinfo)

#now add phenotype info to the covariate file

dfcovarcols <- data.frame(colnames(rev_covar_PC20))
dfcovarcols <- data.frame(dfcovarcols[2:259,])
colnames(dfcovarcols) <- 'genotype'


phenoPC <- read.csv('Detection_of_eQTL/PC_phenotypes.csv',header = T)
sample_rep_batch2 <- read.csv('Detection_of_eQTL/sample_rep_batch2.csv',header = T)
colnames(sample_rep_batch2) <- c('genotype','replicate')


dfcovarcols <- dfcovarcols %>%
  left_join(sample_rep_batch2) %>%
  unite("both",genotype:replicate,sep = "_",remove = F)
dfcovarcols$both <- gsub('_NA','',dfcovarcols$both)
#correction for a few samples
dfcovarcols$both <- gsub('PC_67_4','PC_67_3',dfcovarcols$both)
dfcovarcols$both <- gsub('PC_94_4','PC_94_3',dfcovarcols$both)
dfcovarcols$both <- gsub('PC_15_3_2','PC_15_3',dfcovarcols$both)
dfcovarcols$both <- gsub('PC_16_1_2','PC_16_1',dfcovarcols$both)
dfcovarcols$both <- gsub('PC_23_1_2','PC_23_1',dfcovarcols$both)
dfcovarcols$both <- gsub('PC_26_4_2','PC_26_4',dfcovarcols$both)
dfcovarcols$both <- gsub('PC_35_1_2','PC_35_1',dfcovarcols$both)
dfcovarcols$both <- gsub('PC_38_3_2','PC_38_3',dfcovarcols$both)
dfcovarcols$both <- gsub('PC_46_1_2','PC_46_1',dfcovarcols$both)
dfcovarcols$both <- gsub('PC_50_2_2','PC_50_2',dfcovarcols$both)
dfcovarcols$both <- gsub('PC_53_2_2','PC_53_2',dfcovarcols$both)

phenoPC <- phenoPC %>%
  unite('both',genotype:replicate,sep = "_",remove = F)

joined_covar_phenoPC <- dfcovarcols %>%
  inner_join(phenoPC,by='both') %>%
  dplyr::select(c('genotype.x','Column'))

joinedPCwide <- data.frame(t(joined_covar_phenoPC))
colnames(joinedPCwide) <- joinedPCwide[1,]
joinedPCwide <- joinedPCwide[2,]
joinedPCwide$id <- rownames(joinedPCwide)
rownames(joinedPCwide) <- NULL

#order
joinedPCwide <- joinedPCwide %>% 
  select(id,everything())

rev_covar_PC20_all <- rbind(rev_covar_PC20,joinedPCwide)


#After discussion put NA values in phenotyping data for problematic samples (rep name does not correspond with phenotyping data)
rev_covar_PC20_all[22,c('PC_67','PC_94')] <- NA
write.table(rev_covar_PC20_all,file = 'Detection_of_eQTL/new_prep_data/colum_cov_ready_rev_covar_PC20',quote = F,row.names = F,sep = '\t')


#now filter expression data based on samples in GT file

rev_cpmPC <- read.table('Detection_of_eQTL/new_prep_data/review_cpmf6combatseq_edata',header = T)

revPCallF6HC <- rev_cpmPC %>%
  dplyr::select(matches('PC'))
genesfiltePC <- unique(filteQTLPC$gene.id)
filtePCallF6HC <- revPCallF6HC[which(rownames(revPCallF6HC) %in% genesfiltePC),]

samplesfilteQTLPC <- colnames(corr_rev_both_batchesPC[,2:259])
filtePCallF6HC <- filtePCallF6HC[,which(colnames(filtePCallF6HC) %in% samplesfilteQTLPC)]
filtePCallF6HC <- filtePCallF6HC[,]


filtePCallF6HC_order <- order(colnames(filtePCallF6HC))
filtePCallF6HC <- filtePCallF6HC[,filtePCallF6HC_order]
filtePCallF6HC[filtePCallF6HC==0] <- NA
filtePCallF6HC <- filtePCallF6HC %>%
  dplyr::select(!matches(c('PC_P','PC_C')))
filtePCallF6HC$id <- rownames(filtePCallF6HC)
filtePCallF6HC <- filtePCallF6HC %>%
  select(id,everything())

write.table(filtePCallF6HC,file = 'Detection_of_eQTL/new_prep_data/ready_rev_cpm_combatseqPC',sep = '\t',row.names = F,quote = F)


###########
filteQTLPC <- read.table('Detection_of_eQTL/gatknewway_filteQTLPCC3',header = T)
exgenes <- filteQTLPC[,c('id','gene.id')]


parPCeqtl <- read.table('Detection_of_eQTL/parPCeqtl',header = T)
parPCeqtl[parPCeqtl=="0/0"] <- as.numeric(0)
parPCeqtl[parPCeqtl=="0/1"] <- as.numeric(1)
parPCeqtl[parPCeqtl=="1/1"] <- as.numeric(2)
parPCeqtl[parPCeqtl=="."] <- NA

snpfile_pc <- read.table('Detection_of_eQTL/new_prep_data/ready_rev_gt_PC_both_batches',header = T)
exprfile_pc <- read.table('Detection_of_eQTL/new_prep_data/ready_rev_cpm_combatseqPC',header = T)
covfile_pc <- read.table('Detection_of_eQTL/new_prep_data/colum_cov_ready_rev_covar_PC2_repl',header = T)



snpparPC <- snpfile_pc %>%
  inner_join(parPCeqtl) %>%
  select(id,gene.id,everything())
results <- apply(snpparPC[, 3:260], 2, function(col) {
  ifelse(col == snpparPC[, 261], "P1",
         ifelse(col == snpparPC[, 262], "P2",
                ifelse(col == 1, "HET", NA)))
})

results <- as.data.frame(results)
results[results=="P1"] <- as.numeric(0)
results[results=="HET"] <- as.numeric(1)
results[results=="P2"] <- as.numeric(2)
results <- results %>%
  mutate_all(function(x) as.numeric(as.character(x)))

idgenes <- parPCeqtl[,c('gene.id')]
resultsids <- cbind(idgenes,results)
colnames(resultsids)[1] <- 'gene.id'

uni_resultsids <- resultsids %>%
  group_by(gene.id) %>%
  summarise(across(everything(),~mean(.,na.rm=T)))
uni_resultsids[uni_resultsids<0.5] <- 0
uni_resultsids[uni_resultsids>=0.5 & uni_resultsids <=1] <- 1
uni_resultsids[uni_resultsids<1.5 & uni_resultsids >= 1] <- 1
uni_resultsids[uni_resultsids>=1.5 & uni_resultsids < 2] <- 2


snpupd_pc <- snpfile_pc %>%
  inner_join(exgenes) %>%
  select(id,gene.id,everything()) %>%
  pivot_longer(!id:gene.id,names_to = 'f5line',values_to = 'snp_genotype')

exinugenes <- filteQTLPC %>%
  select('gene.id') %>%
  distinct()
uni_resultsids <- uni_resultsids %>%
  inner_join(exinugenes) %>%
  select(gene.id,everything()) %>%
  pivot_longer(!gene.id,names_to = 'f5line',values_to = 'snp_genotype')


exprupd_pc <- exprfile_pc %>%
  pivot_longer(!id,names_to = 'f5line',values_to = 'cpm_expression')
colnames(exprupd_pc)[1] <- 'gene.id'

tcovfile_pc <- as.data.frame(t(covfile_pc))
colnames(tcovfile_pc) <- tcovfile_pc[1,]
tcovfile_pc <- tcovfile_pc[2:161,]
tcovfile_pc$f5line <- rownames(tcovfile_pc)
rownames(tcovfile_pc) <- NULL


allupd_pc <- uni_resultsids %>%
  inner_join(exprupd_pc) %>%
  inner_join(tcovfile_pc)



allupd_pc$batch <- as.numeric(allupd_pc$batch)
allupd_pc$replicate <- as.numeric(allupd_pc$replicate)
lmm <- lmer(log(cpm_expression) ~ (1 | snp_genotype) + batch + Column + replicate, REML = TRUE,data = allupd_pc)
residuals <- residuals(lmm)
plot(fitted(lmm),residuals)
abline(0,0)
plot(fitted(lmm),residuals,ylim = c(-5,5))

qqnorm(residuals)


model_data <- model.frame(lmm)
allupd_pc$row_id <- seq_len(nrow(allupd_pc))
used_rows <- rownames(model_data)

# Create a vector of NA values matching the original data's length
allupd_pc$residuals <- NA

residdf <- as.data.frame(residuals(lmm))

# Fill in the residuals for rows used in the model
allupd_pc$residuals[as.integer(used_rows)] <- residdf[, 1]

# write.table(allupd_pc,'Detection_of_eQTL/new_prep_data/post_lmmPC_singlegene',quote = F,row.names = F,sep = '\t')
########
# Now split again for the eQTL analysis into Genotype and Expression file with residuals as values
nonNApostlmm <- allupd_pc[!is.na(allupd_pc$residuals),]

### for single gene approach

expr_postlmmPC <- nonNApostlmm[,c('gene.id','f5line','residuals')]
expr_postlmmPC <- expr_postlmmPC %>%
  pivot_wider(names_from = 'f5line',values_from = 'residuals')
write.table(expr_postlmmPC,'Detection_of_eQTL/new_prep_data/postlmm_expr_PCsinglegene',quote = F,row.names = F,sep = '\t')

snp_postlmmPC <- nonNApostlmm %>%
  select(gene.id) %>%
  distinct() %>%
  inner_join(snpupd_pc,relationship = 'many-to-many') %>%
  pivot_wider(names_from = 'f5line',values_from = 'snp_genotype') %>%
  select(!gene.id)

namevec <- unlist(nonNApostlmm[,c('f5line')])
namevec <- c('id',namevec)
namevec <- unlist(colnames(expr_postlmmPC[2:159]))
namevec <- c('id',namevec)

snp_postlmmPC <- snp_postlmmPC[,names(snp_postlmmPC) %in% namevec]
snp_postlmmPC <- snp_postlmmPC[,namevec]


write.table(snp_postlmmPC,'Detection_of_eQTL/new_prep_data/postlmm_genotype_PCsinglegene',quote = F,row.names = F,sep = '\t')



#Now rerun eQTL with the new files
library(MatrixEQTL)
## Settings
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR;

# Genotype file name
SNP_file_name = paste("Detection_of_eQTL/new_prep_data/postlmm_genotype_PCsinglegene", sep="");

# Gene expression file name
expression_file_name = paste("Detection_of_eQTL/new_prep_data/postlmm_expr_PCsinglegene", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = character();

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-5;
pvOutputThreshold_tra = 1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

snpspos = read.table("Detection_of_eQTL/newway_snp_locPCC3", header = TRUE, stringsAsFactors = FALSE);
genepos = read.table("Detection_of_eQTL/newway_gene_locPCC3", header = TRUE, stringsAsFactors = FALSE);

#reduce because the nmber of SNPs and genes is reduced after the LMM
snpsposupd <- snpspos[snpspos$id %in% unlist(snp_postlmmPC[,c('id')]),]
geneposupd <- genepos[genepos$gene.id %in% unlist(expr_postlmmPC[,c('gene.id')]),]

## Run the analysis

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpsposupd,
  genepos = geneposupd,
  cisDist = cisDist,
  pvalue.hist = 'qqplot',
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

#### Results:
plot(me)



mecis <- me$cis$eqtls
metrans <- me$trans$eqtls

### Chr bias?
geneallloc <- geneposupd
colnames(geneallloc)[1] <- c('gene')
cischrbias <- geneallloc %>%
  inner_join(mecis,relationship = 'many-to-many')
ggplot(cischrbias,aes(chr)) + geom_bar()

transchrbias <- geneallloc %>%
  inner_join(metrans,relationship = 'many-to-many')
ggplot(transchrbias,aes(chr)) + geom_bar()

filteQTLPC <- read.table('Detection_of_eQTL/gatknewway_filteQTLPCC3',header = T)

basePCsnpinfo <- filteQTLPC[,1:6]
colnames(basePCsnpinfo) <- c('snps','chr','pos','gene.id','start','end')
basePCsnpinfo$snps <- as.character(basePCsnpinfo$snps)
mergePCcis <- basePCsnpinfo %>%
  inner_join(mecis,relationship = "many-to-many")

mergePCctrans <- basePCsnpinfo %>%
  inner_join(metrans,relationship = "many-to-many")


##
colnames(mergePCcis)[4] <- 'snp_gene_origin'
added_pheno_PCcis <-mergePCcis

colnames(mergePCctrans)[4] <- 'snp_gene_origin'
added_pheno_PCtrans <-mergePCctrans
added_pheno_PCtrans$gene_affected_chr <- substr(added_pheno_PCtrans$gene,start=8,stop = 9)
added_pheno_PCtrans_filt_diff_chr <- added_pheno_PCtrans[added_pheno_PCtrans$chr!=added_pheno_PCtrans$gene_affected_chr,]
added_pheno_PCtrans_filt_same_chr <- added_pheno_PCtrans[added_pheno_PCtrans$chr==added_pheno_PCtrans$gene_affected_chr,]


# Now do the filtering for FDR

added_pheno_PCcis <- added_pheno_PCcis[added_pheno_PCcis$FDR<0.001,]
added_pheno_PCtrans_filt_diff_chr <- added_pheno_PCtrans_filt_diff_chr[added_pheno_PCtrans_filt_diff_chr$FDR<0.001,]
added_pheno_PCtrans_filt_same_chr <- added_pheno_PCtrans_filt_same_chr[added_pheno_PCtrans_filt_same_chr$FDR<0.001,]

#
#
#
#
# Now do the UPset plots to compare results with results in the preprint

# Output new results
write.table(added_pheno_PCcis,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PCcis_FDR001',row.names = F,quote = F,sep = '\t')
write.table(added_pheno_PCtrans_filt_diff_chr,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PCtrans_filt_diff_chr_FDR001',row.names = F,quote = F,sep = '\t')
write.table(added_pheno_PCtrans_filt_same_chr,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PCtrans_filt_same_chr_FDR001',row.names = F,quote = F,sep = '\t')



# Continue doing other stats


added_pheno_PCtrans_filt_diff_chr <- added_pheno_PCtrans_filt_diff_chr[,1:7]
added_pheno_PCtrans_filt_same_chr <- added_pheno_PCtrans_filt_same_chr[,1:7]
added_pheno_PCcis <- added_pheno_PCcis[,1:7]


# SNP inflation is caused by multiple SNPs having association with the same gene - we only want to see one shared example
# So we will get rid of the SNP id and its position information

added_pheno_PCcis <- added_pheno_PCcis %>%
  dplyr::select(!c('snps','pos'))
added_pheno_PCcis <- unique(added_pheno_PCcis)
added_pheno_PCtrans_filt_diff_chr <- added_pheno_PCtrans_filt_diff_chr %>%
  dplyr::select(!c('snps','pos'))
added_pheno_PCtrans_filt_diff_chr <- unique(added_pheno_PCtrans_filt_diff_chr)
added_pheno_PCtrans_filt_same_chr <- added_pheno_PCtrans_filt_same_chr %>%
  dplyr::select(!c('snps','pos'))
added_pheno_PCtrans_filt_same_chr <- unique(added_pheno_PCtrans_filt_same_chr)

# The numbers reduced significantly and now we can make some stats

# First lets check the numbers for each genome

added_pheno_PCcis <- added_pheno_PCcis %>%
  mutate(sub_gene_origin = substr(snp_gene_origin,9,9)) %>%
  mutate(sub_gene_regulated = substr(gene,9,9))
added_pheno_PCtrans_filt_diff_chr <- added_pheno_PCtrans_filt_diff_chr %>%
  mutate(sub_gene_origin = substr(snp_gene_origin,9,9)) %>%
  mutate(sub_gene_regulated = substr(gene,9,9))
added_pheno_PCtrans_filt_same_chr <- added_pheno_PCtrans_filt_same_chr %>%
  mutate(sub_gene_origin = substr(snp_gene_origin,9,9)) %>%
  mutate(sub_gene_regulated = substr(gene,9,9))

write.table(added_pheno_PCcis,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PCcis_deflated',row.names = F,sep = '\t',quote = F)
write.table(added_pheno_PCtrans_filt_diff_chr,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PCtrans_filt_diff_chr_deflated',row.names = F,sep = '\t',quote = F)
write.table(added_pheno_PCtrans_filt_same_chr,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PCtrans_filt_same_chr_deflated',row.names = F,sep = '\t',quote = F)


# Make graphs summarizing the eQTL results using upset plots
library(ComplexUpset)

PCcis <- read.table('Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PCcis_deflated',header = T)
PCtrans_filt_diff_chr <- read.table('Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PCtrans_filt_diff_chr_deflated',header=T)
PCtrans_filt_same_chr <- read.table('Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PCtrans_filt_same_chr_deflated',header = T)

PCcis <- added_pheno_PCcis
PCtrans_filt_diff_chr <- added_pheno_PCtrans_filt_diff_chr
PCtrans_filt_same_chr <- added_pheno_PCtrans_filt_same_chr


PCtransdiff_upset <- PCtrans_filt_diff_chr
PCtransdiff_upset <- PCtransdiff_upset %>%
  mutate(snp_A = as.numeric(sub_gene_origin == 'A'),
         snp_B = as.numeric(sub_gene_origin == 'B'),
         snp_D = as.numeric(sub_gene_origin == 'D'),
         gene_A = as.numeric(sub_gene_regulated == 'A'),
         gene_B = as.numeric(sub_gene_regulated == 'B'),
         gene_D = as.numeric(sub_gene_regulated == 'D'),
  )

PCtranssame_upset <- PCtrans_filt_same_chr
PCtranssame_upset <- PCtranssame_upset %>%
  mutate(snp_A = as.numeric(sub_gene_origin == 'A'),
         snp_B = as.numeric(sub_gene_origin == 'B'),
         snp_D = as.numeric(sub_gene_origin == 'D'),
         gene_A = as.numeric(sub_gene_regulated == 'A'),
         gene_B = as.numeric(sub_gene_regulated == 'B'),
         gene_D = as.numeric(sub_gene_regulated == 'D'),
  )
PCcis_upset <- PCcis
PCcis_upset <- PCcis_upset %>%
  mutate(snp_A = as.numeric(sub_gene_origin == 'A'),
         snp_B = as.numeric(sub_gene_origin == 'B'),
         snp_D = as.numeric(sub_gene_origin == 'D'),
         gene_A = as.numeric(sub_gene_regulated == 'A'),
         gene_B = as.numeric(sub_gene_regulated == 'B'),
         gene_D = as.numeric(sub_gene_regulated == 'D'),
  )


cispc <- upset(PCcis_upset,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
               sort_intersections=F,
               intersections=list(c('snp_A','gene_A'),c('snp_B','gene_B'),c('snp_D','gene_D')),
               stripes = c(alpha("#FF950E",0.9),alpha("#4B1F6F",0.9),alpha("#579D1C",0.9)),
               matrix=intersection_matrix(geom = geom_point(size=7),segment = geom_segment(lwd=4)),
               base_annotations = list('size'=(intersection_size(width=0.3,text = list(size=9),bar_number_threshold = 1) + 
                                                 theme_classic() +
                                                 theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                                       axis.title = element_text(size = 23),
                                                       title = element_text(size = 18)) +
                                                 ylim(0,10000) +
                                                 ylab('Number of associations') + xlab("") +
                                                 ggtitle('PC cis-eQTL associations'))),
               width_ratio = 0.2,queries = list(
                 upset_query(intersect=c('snp_A','gene_A'),fill='#579D1C'),
                 upset_query(intersect=c('snp_B','gene_B'),fill='#4B1F6F'),
                 upset_query(intersect=c('snp_D','gene_D'),fill='#FF950E'),
                 upset_query(set='snp_A', fill='black'),
                 upset_query(set='snp_B', fill='black'),
                 upset_query(set='snp_D', fill='black'),
                 upset_query(set='gene_A', fill='black'),
                 upset_query(set='gene_B', fill='black'),
                 upset_query(set='gene_D', fill='black')
                 
               )) + xlab("") + theme(axis.text = element_text(size = 18))
cispc

transspc <- upset(PCtranssame_upset,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
                  sort_intersections=F,
                  intersections=list(c('snp_A','gene_A'),c('snp_B','gene_B'),c('snp_D','gene_D')),
                  stripes = c(alpha("#FF950E",0.9),alpha("#4B1F6F",0.9),alpha("#579D1C",0.9)),
                  matrix=intersection_matrix(geom = geom_point(size=7),segment = geom_segment(lwd=4)),
                  base_annotations = list(
                    'size'=(intersection_size(width=0.3,text = list(size=9)) + 
                              theme_classic() +
                              theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                    axis.title = element_text(size = 16),
                                    title = element_text(size = 18)) +
                              coord_cartesian(ylim = c(0,10000)) +
                              ylab('Number of associations') + xlab(""))
                  ),
                  width_ratio = 0.2,queries = list(
                    upset_query(intersect=c('snp_A','gene_A'),fill='#579D1C'),
                    upset_query(intersect=c('snp_B','gene_B'),fill='#4B1F6F'),
                    upset_query(intersect=c('snp_D','gene_D'),fill='#FF950E'),
                    upset_query(set='snp_A', fill='black'),
                    upset_query(set='snp_B', fill='black'),
                    upset_query(set='snp_D', fill='black'),
                    upset_query(set='gene_A', fill='black'),
                    upset_query(set='gene_B', fill='black'),
                    upset_query(set='gene_D', fill='black')
                    
                  )) + xlab("") + theme(axis.text = element_text(size = 18))
transspc

transspctop <- upset(PCtranssame_upset,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
                     sort_intersections=F,
                     intersections=list(c('snp_A','gene_A'),c('snp_B','gene_B'),c('snp_D','gene_D')),
                     stripes = c(alpha("#FF950E",0.9),alpha("#4B1F6F",0.9),alpha("#579D1C",0.9)),
                     matrix=intersection_matrix(geom = geom_point(size=7),segment = geom_segment(lwd=4)),
                     base_annotations = list(
                       'topsize'=(intersection_size(width=0.3,text = list(size=9),bar_number_threshold = 1) + 
                                    theme_classic() +
                                    theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                          axis.title = element_text(size = 23),
                                          title = element_text(size = 18)) +
                                    coord_cartesian(ylim = c(20000,55000)) +
                                    ylab('') + xlab("") +
                                    ggtitle('PC trans-eQTL same chromosome associations'))
                     ),
                     width_ratio = 0.2,queries = list(
                       upset_query(intersect=c('snp_A','gene_A'),fill='#579D1C'),
                       upset_query(intersect=c('snp_B','gene_B'),fill='#4B1F6F'),
                       upset_query(intersect=c('snp_D','gene_D'),fill='#FF950E'),
                       upset_query(set='snp_A', fill='black'),
                       upset_query(set='snp_B', fill='black'),
                       upset_query(set='snp_D', fill='black'),
                       upset_query(set='gene_A', fill='black'),
                       upset_query(set='gene_B', fill='black'),
                       upset_query(set='gene_D', fill='black')
                       
                     )) + xlab("") + theme(axis.text = element_text(size = 18))
transspctop

transdpc <- upset(PCtransdiff_upset,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
                  sort_intersections=F,
                  intersections=list(c('snp_A','gene_A'),c('snp_A','gene_B'),c('snp_A','gene_D'),
                                     c('snp_B','gene_A'),c('snp_B','gene_B'),c('snp_B','gene_D'),
                                     c('snp_D','gene_A'),c('snp_D','gene_B'),c('snp_D','gene_D')),
                  stripes = c(alpha("#FF950E",0.9),alpha("#4B1F6F",0.9),alpha("#579D1C",0.9)),
                  matrix=intersection_matrix(geom = geom_point(size=7),segment = geom_segment(lwd=4)),
                  base_annotations = list('size'=(intersection_size(width=0.9,text = list(size=9),bar_number_threshold = 1) + 
                                                    theme_classic() +
                                                    theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                                          axis.title = element_text(size = 23),
                                                          title = element_text(size = 18)) +
                                                    ylab('Number of associations') + xlab("") +
                                                    ylim(0,10000) +
                                                    ggtitle('PC trans-eQTL different chromosome associations'))),
                  width_ratio = 0.2,queries = list(
                    upset_query(intersect=c('snp_A','gene_A'),fill='#579D1C'),
                    upset_query(intersect=c('snp_B','gene_B'),fill='#4B1F6F'),
                    upset_query(intersect=c('snp_D','gene_D'),fill='#FF950E'),
                    upset_query(intersect=c('snp_A','gene_D'),fill='black'),
                    upset_query(intersect=c('snp_B','gene_D'),fill='black'),
                    upset_query(intersect=c('snp_A','gene_B'),fill='black'),
                    upset_query(intersect=c('snp_B','gene_A'),fill='black'),
                    upset_query(intersect=c('snp_D','gene_A'),fill='black'),
                    upset_query(intersect=c('snp_D','gene_B'),fill='black')
                    
                  )) + xlab("") + theme(axis.text = element_text(size = 18))
transdpc

ggsave(plot = transdpc,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transdpc_upset_associations.png',device = 'png',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = transspc,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transspc_upset_associations.png',device = 'png',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = cispc,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/cispc_upset_associations.png',device = 'png',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = transspctop,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transspctop_upset_associations.png',device = 'png',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = transdpc,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transdpc_upset_associations.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = transspc,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transspc_upset_associations.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = cispc,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/cispc_upset_associations.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = transspctop,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transspctop_upset_associations.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)







# PxW nem model eQTL


filteQTLPW <- read.table('Detection_of_eQTL/gatknewway_filteQTLPW_batch1',header = T)

rev_gatknewwayPW <- read.table('Detection_of_eQTL/gatknewway_gt_PW_batch1',header = T)

rev_pcallbatchdf <- read.table('Detection_of_eQTL/new_prep_data/f6allbatchdf',header = T) 
rev_pcallbatchdf <- rev_pcallbatchdf[rev_pcallbatchdf$batch==1 & rev_pcallbatchdf$cross==2,]
rev_pcallbatchdf_cols <- unlist(rownames(rev_pcallbatchdf))
rev_pcallbatchdf_cols <- c('id',rev_pcallbatchdf_cols)

rev_gatknewwayPW <- rev_gatknewwayPW[,names(rev_gatknewwayPW) %in% rev_pcallbatchdf_cols]
write.table(rev_gatknewwayPW,'Detection_of_eQTL/new_prep_data/rev_gatknewwayPW_excel',quote = F,row.names = F,sep = '\t')


# Load in new GT data and create for covariate file
corr_rev_gatknewwayPW <- read.table('Detection_of_eQTL/new_prep_data/correct_rev_gatknewwayPW_excel',header = T)

#now add batch info to the covariate file

batchinfo <- read.table('Detection_of_eQTL/new_prep_data/f6allbatchdf',header = T)
batchinfo <- batchinfo[batchinfo$batch==1 & batchinfo$cross==2,]
batchinfo <- as.data.frame(t(batchinfo))
batchinfo$id <- rownames(batchinfo)
batchinfo <- batchinfo %>%
  select(id,everything())
rownames(batchinfo) <- NULL
batchinfo <- batchinfo %>%
  filter(!row_number() %in% c(2))


rev_covarPW <- batchinfo[,1:151]

#now add phenotype info to the covariate file

dfcovarcols <- data.frame(colnames(rev_covarPW))
dfcovarcols <- data.frame(dfcovarcols[2:157,])
colnames(dfcovarcols) <- 'samples'


phenoPW <- read.csv('Detection_of_eQTL/PW_phenotypes.csv',header = T)

phenoPW$samples <- gsub('-','_',phenoPW$samples)

joined_covar_phenoPW <- dfcovarcols %>%
  inner_join(phenoPW,by='samples') %>%
  dplyr::select(c('samples','Column'))

joinedPWwide <- data.frame(t(joined_covar_phenoPW))
colnames(joinedPWwide) <- joinedPWwide[1,]
joinedPWwide <- joinedPWwide[2,]
joinedPWwide$id <- rownames(joinedPWwide)
rownames(joinedPWwide) <- NULL


#order
joinedPWwide <- joinedPWwide %>% 
  select(id,everything())

rev_covar_PW_all <- rbind(rev_covarPW,joinedPWwide)
write.table(rev_covar_PW_all,file = 'Detection_of_eQTL/new_prep_data/colum_cov_ready_rev_covar_PW',quote = F,row.names = F,sep = '\t')
rev_covar_PW_all <- read.table('Detection_of_eQTL/new_prep_data/colum_cov_ready_rev_covar_PW',header = T)

#now filter expression data based on samples in GT file

rev_cpmPW <- read.table('Detection_of_eQTL/new_prep_data/review_cpmf6combatseq_edata',header = T)

revPWallF6HC <- rev_cpmPW %>%
  dplyr::select(matches('PW'))
genesfiltePW <- unique(filteQTLPW$gene.id)
filtePWallF6HC <- revPWallF6HC[which(rownames(revPWallF6HC) %in% genesfiltePW),]

samplesfilteQTLPW <- colnames(rev_covar_PW_all[,2:151])
filtePWallF6HC <- filtePWallF6HC[,which(colnames(filtePWallF6HC) %in% samplesfilteQTLPW)]

filtePWallF6HC[filtePWallF6HC==0] <- NA
filtePWallF6HC$id <- rownames(filtePWallF6HC)
filtePWallF6HC <- filtePWallF6HC %>%
  select(id,everything())

write.table(filtePWallF6HC,file = 'Detection_of_eQTL/new_prep_data/ready_rev_cpm_combatseqPW',sep = '\t',row.names = F,quote = F)


###########
filteQTLPW <- read.table('Detection_of_eQTL/gatknewway_filteQTLPW_batch1',header = T)
exgenes <- filteQTLPW[,c('id','gene.id')]


parPWeqtl <- read.table('Detection_of_eQTL/parPWeqtl',header = T)
parPWeqtl[parPWeqtl=="0/0"] <- as.numeric(0)
parPWeqtl[parPWeqtl=="0/1"] <- as.numeric(1)
parPWeqtl[parPWeqtl=="1/1"] <- as.numeric(2)
parPWeqtl[parPWeqtl=="."] <- NA

snpfile_pw <- read.table('Detection_of_eQTL/new_prep_data/correct_rev_gatknewwayPW_excel',header = T)

exprfile_pw <- read.table('Detection_of_eQTL/new_prep_data/ready_rev_cpm_combatseqPW',header = T)

covfile_pw <- read.table('Detection_of_eQTL/new_prep_data/colum_cov_ready_rev_covar_PW_repl',header = T)


snpparPW <- snpfile_pw %>%
  inner_join(parPWeqtl) %>%
  select(id,gene.id,everything())
results <- apply(snpparPW[, 3:152], 2, function(col) {
  ifelse(col == snpparPW[, 153], "P1",
         ifelse(col == snpparPW[, 154], "P2",
                ifelse(col == 1, "HET", NA)))
})



results <- as.data.frame(results)
results[results=="P1"] <- as.numeric(0)
results[results=="HET"] <- as.numeric(1)
results[results=="P2"] <- as.numeric(2)
results <- results %>%
  mutate_all(function(x) as.numeric(as.character(x)))

idgenes <- parPWeqtl[,c('gene.id')]
resultsids <- cbind(idgenes,results)
colnames(resultsids)[1] <- 'gene.id'

uni_resultsids <- resultsids %>%
  group_by(gene.id) %>%
  summarise(across(everything(),~mean(.,na.rm=T)))
uni_resultsids[uni_resultsids<0.5] <- 0
uni_resultsids[uni_resultsids>=0.5 & uni_resultsids <=1] <- 1
uni_resultsids[uni_resultsids<1.5 & uni_resultsids >= 1] <- 1
uni_resultsids[uni_resultsids>=1.5 & uni_resultsids < 2] <- 2


snpupd_pw <- snpfile_pw %>%
  inner_join(exgenes) %>%
  select(id,gene.id,everything()) %>%
  pivot_longer(!id:gene.id,names_to = 'f5line',values_to = 'snp_genotype')

exinugenes <- filteQTLPW %>%
  select('gene.id') %>%
  distinct()
uni_resultsids <- uni_resultsids %>%
  inner_join(exinugenes) %>%
  select(gene.id,everything()) %>%
  pivot_longer(!gene.id,names_to = 'f5line',values_to = 'snp_genotype')


exprupd_pw <- exprfile_pw %>%
  pivot_longer(!id,names_to = 'f5line',values_to = 'cpm_expression')
colnames(exprupd_pw)[1] <- 'gene.id'

tcovfile_pw <- as.data.frame(t(covfile_pw))
colnames(tcovfile_pw) <- tcovfile_pw[1,]
tcovfile_pw <- tcovfile_pw[2:151,]
tcovfile_pw$f5line <- rownames(tcovfile_pw)
rownames(tcovfile_pw) <- NULL


allupd_pw <- uni_resultsids %>%
  inner_join(exprupd_pw) %>%
  inner_join(tcovfile_pw)


library(lme4)
allupd_pw$batch <- as.numeric(allupd_pw$batch)
allupd_pw$replicate <- as.numeric(allupd_pw$replicate)
lmm <- lmer(log(cpm_expression) ~ (1 | snp_genotype) + batch + Column + replicate, REML = TRUE,data = allupd_pw)
residuals <- residuals(lmm)
plot(fitted(lmm),residuals)
abline(0,0)
plot(fitted(lmm),residuals,ylim = c(-5,5))

qqnorm(residuals)


model_data <- model.frame(lmm)
allupd_pw$row_id <- seq_len(nrow(allupd_pw))
used_rows <- rownames(model_data)

# Create a vector of NA values matching the original data's length
allupd_pw$residuals <- NA

residdf <- as.data.frame(residuals(lmm))

# Fill in the residuals for rows used in the model
allupd_pw$residuals[as.integer(used_rows)] <- residdf[, 1]

write.table(allupd_pw,'Detection_of_eQTL/new_prep_data/post_lmmPW_singlegene',quote = F,row.names = F,sep = '\t')
########
# Now split again for the eQTL analysis into Genotype and Expression file with residuals as values
nonNApostlmm <- allupd_pw[!is.na(allupd_pw$residuals),]

### for single gene approach

expr_postlmmPW <- nonNApostlmm[,c('gene.id','f5line','residuals')]
expr_postlmmPW <- expr_postlmmPW %>%
  pivot_wider(names_from = 'f5line',values_from = 'residuals')
write.table(expr_postlmmPW,'Detection_of_eQTL/new_prep_data/postlmm_expr_PWsinglegene',quote = F,row.names = F,sep = '\t')

snp_postlmmPW <- nonNApostlmm %>%
  select(gene.id) %>%
  distinct() %>%
  inner_join(snpupd_pw,relationship = 'many-to-many') %>%
  pivot_wider(names_from = 'f5line',values_from = 'snp_genotype') %>%
  select(!gene.id)

namevec <- unlist(nonNApostlmm[,c('f5line')])
namevec <- c('id',namevec)
namevec <- unlist(colnames(expr_postlmmPW[2:151]))
namevec <- c('id',namevec)

snp_postlmmPW <- snp_postlmmPW[,names(snp_postlmmPW) %in% namevec]
snp_postlmmPW <- snp_postlmmPW[,namevec]


write.table(snp_postlmmPW,'Detection_of_eQTL/new_prep_data/postlmm_genotype_PWsinglegene',quote = F,row.names = F,sep = '\t')


#Now rerun eQTL with the new files
library(MatrixEQTL)
## Settings
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR;

# Genotype file name
SNP_file_name = paste("Detection_of_eQTL/new_prep_data/postlmm_genotype_PWsinglegene", sep="");

# Gene expression file name
expression_file_name = paste("Detection_of_eQTL/new_prep_data/postlmm_expr_PWsinglegene", sep="");

# Covariates file name
# Set to character() for no covariates
covariates_file_name = character();

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1e-5;
pvOutputThreshold_tra = 1e-5;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

snpspos = read.table("Detection_of_eQTL/gatknewway_snp_locPW_batch1", header = TRUE, stringsAsFactors = FALSE);
genepos = read.table("Detection_of_eQTL/gatknewway_gene_locPW_batch1", header = TRUE, stringsAsFactors = FALSE);

#reduce because the nmber of SNPs and genes is reduced after the LMM
snpsposupd <- snpspos[snpspos$id %in% unlist(snp_postlmmPW[,c('id')]),]
geneposupd <- genepos[genepos$gene.id %in% unlist(expr_postlmmPW[,c('gene.id')]),]

## Run the analysis

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name     = output_file_name_tra,
  pvOutputThreshold     = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpsposupd,
  genepos = geneposupd,
  cisDist = cisDist,
  pvalue.hist = 'qqplot',
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_tra);
unlink(output_file_name_cis);

#### Results:
plot(me)



mecis <- me$cis$eqtls
metrans <- me$trans$eqtls

### Chr bias?
geneallloc <- geneposupd
colnames(geneallloc)[1] <- c('gene')
cischrbias <- geneallloc %>%
  inner_join(mecis,relationship = 'many-to-many')
ggplot(cischrbias,aes(chr)) + geom_bar()

transchrbias <- geneallloc %>%
  inner_join(metrans,relationship = 'many-to-many')
ggplot(transchrbias,aes(chr)) + geom_bar()

filteQTLPW <- read.table('Detection_of_eQTL/gatknewway_filteQTLPW_batch1',header = T)

basePWsnpinfo <- filteQTLPW[,1:6]
colnames(basePWsnpinfo) <- c('snps','chr','pos','gene.id','start','end')
basePWsnpinfo$snps <- as.character(basePWsnpinfo$snps)
mergePWcis <- basePWsnpinfo %>%
  inner_join(mecis,relationship = "many-to-many")

mergePWctrans <- basePWsnpinfo %>%
  inner_join(metrans,relationship = "many-to-many")


##
colnames(mergePWcis)[4] <- 'snp_gene_origin'
added_pheno_PWcis <-mergePWcis

colnames(mergePWctrans)[4] <- 'snp_gene_origin'
added_pheno_PWtrans <-mergePWctrans
added_pheno_PWtrans$gene_affected_chr <- substr(added_pheno_PWtrans$gene,start=8,stop = 9)
added_pheno_PWtrans_filt_diff_chr <- added_pheno_PWtrans[added_pheno_PWtrans$chr!=added_pheno_PWtrans$gene_affected_chr,]
added_pheno_PWtrans_filt_same_chr <- added_pheno_PWtrans[added_pheno_PWtrans$chr==added_pheno_PWtrans$gene_affected_chr,]


# Now do the filtering for FDR

added_pheno_PWcis <- added_pheno_PWcis[added_pheno_PWcis$FDR<0.001,]
added_pheno_PWtrans_filt_diff_chr <- added_pheno_PWtrans_filt_diff_chr[added_pheno_PWtrans_filt_diff_chr$FDR<0.001,]
added_pheno_PWtrans_filt_same_chr <- added_pheno_PWtrans_filt_same_chr[added_pheno_PWtrans_filt_same_chr$FDR<0.001,]

#
#
#
#
# Now do the UPset plots to compare results with results in the preprint

# Output new results
write.table(added_pheno_PWcis,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PWcis_FDR001',row.names = F,quote = F,sep = '\t')
write.table(added_pheno_PWtrans_filt_diff_chr,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PWtrans_filt_diff_chr_FDR001',row.names = F,quote = F,sep = '\t')
write.table(added_pheno_PWtrans_filt_same_chr,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PWtrans_filt_same_chr_FDR001',row.names = F,quote = F,sep = '\t')



# Continue doing other stats


added_pheno_PWtrans_filt_diff_chr <- added_pheno_PWtrans_filt_diff_chr[,1:7]
added_pheno_PWtrans_filt_same_chr <- added_pheno_PWtrans_filt_same_chr[,1:7]
added_pheno_PWcis <- added_pheno_PWcis[,1:7]


# SNP inflation is caused by multiple SNPs having association with the same gene - we only want to see one shared example
# So we will get rid of the SNP id and its position information

added_pheno_PWcis <- added_pheno_PWcis %>%
  dplyr::select(!c('snps','pos'))
added_pheno_PWcis <- unique(added_pheno_PWcis)
added_pheno_PWtrans_filt_diff_chr <- added_pheno_PWtrans_filt_diff_chr %>%
  dplyr::select(!c('snps','pos'))
added_pheno_PWtrans_filt_diff_chr <- unique(added_pheno_PWtrans_filt_diff_chr)
added_pheno_PWtrans_filt_same_chr <- added_pheno_PWtrans_filt_same_chr %>%
  dplyr::select(!c('snps','pos'))
added_pheno_PWtrans_filt_same_chr <- unique(added_pheno_PWtrans_filt_same_chr)

# The numbers reduced significantly and now we can make some stats

# First lets check the numbers for each genome

added_pheno_PWcis <- added_pheno_PWcis %>%
  mutate(sub_gene_origin = substr(snp_gene_origin,9,9)) %>%
  mutate(sub_gene_regulated = substr(gene,9,9))
added_pheno_PWtrans_filt_diff_chr <- added_pheno_PWtrans_filt_diff_chr %>%
  mutate(sub_gene_origin = substr(snp_gene_origin,9,9)) %>%
  mutate(sub_gene_regulated = substr(gene,9,9))
added_pheno_PWtrans_filt_same_chr <- added_pheno_PWtrans_filt_same_chr %>%
  mutate(sub_gene_origin = substr(snp_gene_origin,9,9)) %>%
  mutate(sub_gene_regulated = substr(gene,9,9))

write.table(added_pheno_PWcis,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PWcis_deflated',row.names = F,sep = '\t',quote = F)
write.table(added_pheno_PWtrans_filt_diff_chr,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PWtrans_filt_diff_chr_deflated',row.names = F,sep = '\t',quote = F)
write.table(added_pheno_PWtrans_filt_same_chr,'Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PWtrans_filt_same_chr_deflated',row.names = F,sep = '\t',quote = F)


# Make graphs summarizing the eQTL results using upset plots
library(ComplexUpset)

PWcis <- read.table('Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PWcis_deflated',header = T)
PWtrans_filt_diff_chr <- read.table('Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PWtrans_filt_diff_chr_deflated',header=T)
PWtrans_filt_same_chr <- read.table('Detection_of_eQTL/new_prep_data/post_lmm_added_pheno_PWtrans_filt_same_chr_deflated',header = T)

PWcis <- added_pheno_PWcis
PWtrans_filt_diff_chr <- added_pheno_PWtrans_filt_diff_chr
PWtrans_filt_same_chr <- added_pheno_PWtrans_filt_same_chr


PWtransdiff_upset <- PWtrans_filt_diff_chr
PWtransdiff_upset <- PWtransdiff_upset %>%
  mutate(snp_A = as.numeric(sub_gene_origin == 'A'),
         snp_B = as.numeric(sub_gene_origin == 'B'),
         snp_D = as.numeric(sub_gene_origin == 'D'),
         gene_A = as.numeric(sub_gene_regulated == 'A'),
         gene_B = as.numeric(sub_gene_regulated == 'B'),
         gene_D = as.numeric(sub_gene_regulated == 'D'),
  )

PWtranssame_upset <- PWtrans_filt_same_chr
PWtranssame_upset <- PWtranssame_upset %>%
  mutate(snp_A = as.numeric(sub_gene_origin == 'A'),
         snp_B = as.numeric(sub_gene_origin == 'B'),
         snp_D = as.numeric(sub_gene_origin == 'D'),
         gene_A = as.numeric(sub_gene_regulated == 'A'),
         gene_B = as.numeric(sub_gene_regulated == 'B'),
         gene_D = as.numeric(sub_gene_regulated == 'D'),
  )
PWcis_upset <- PWcis
PWcis_upset <- PWcis_upset %>%
  mutate(snp_A = as.numeric(sub_gene_origin == 'A'),
         snp_B = as.numeric(sub_gene_origin == 'B'),
         snp_D = as.numeric(sub_gene_origin == 'D'),
         gene_A = as.numeric(sub_gene_regulated == 'A'),
         gene_B = as.numeric(sub_gene_regulated == 'B'),
         gene_D = as.numeric(sub_gene_regulated == 'D'),
  )


cispw <- upset(PWcis_upset,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
               sort_intersections=F,
               intersections=list(c('snp_A','gene_A'),c('snp_B','gene_B'),c('snp_D','gene_D')),
               stripes = c(alpha("#FF950E",0.9),alpha("#4B1F6F",0.9),alpha("#579D1C",0.9)),
               matrix=intersection_matrix(geom = geom_point(size=7),segment = geom_segment(lwd=4)),
               base_annotations = list('size'=(intersection_size(width=0.3,text = list(size=9),bar_number_threshold = 1) + 
                                                 theme_classic() +
                                                 theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                                       axis.title = element_text(size = 23),
                                                       title = element_text(size = 18)) +
                                                 ylim(0,10000) +
                                                 ylab('Number of associations') + xlab("") +
                                                 ggtitle('PW cis-eQTL associations'))),
               width_ratio = 0.2,queries = list(
                 upset_query(intersect=c('snp_A','gene_A'),fill='#579D1C'),
                 upset_query(intersect=c('snp_B','gene_B'),fill='#4B1F6F'),
                 upset_query(intersect=c('snp_D','gene_D'),fill='#FF950E'),
                 upset_query(set='snp_A', fill='black'),
                 upset_query(set='snp_B', fill='black'),
                 upset_query(set='snp_D', fill='black'),
                 upset_query(set='gene_A', fill='black'),
                 upset_query(set='gene_B', fill='black'),
                 upset_query(set='gene_D', fill='black')
                 
               )) + xlab("") + theme(axis.text = element_text(size = 18))
cispw

transspw <- upset(PWtranssame_upset,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
                  sort_intersections=F,
                  intersections=list(c('snp_A','gene_A'),c('snp_B','gene_B'),c('snp_D','gene_D')),
                  stripes = c(alpha("#FF950E",0.9),alpha("#4B1F6F",0.9),alpha("#579D1C",0.9)),
                  matrix=intersection_matrix(geom = geom_point(size=7),segment = geom_segment(lwd=4)),
                  base_annotations = list(
                    'size'=(intersection_size(width=0.3,text = list(size=9)) + 
                              theme_classic() +
                              theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                    axis.title = element_text(size = 16),
                                    title = element_text(size = 18)) +
                              coord_cartesian(ylim = c(0,15000)) +
                              ylab('Number of associations') + xlab(""))
                  ),
                  width_ratio = 0.2,queries = list(
                    upset_query(intersect=c('snp_A','gene_A'),fill='#579D1C'),
                    upset_query(intersect=c('snp_B','gene_B'),fill='#4B1F6F'),
                    upset_query(intersect=c('snp_D','gene_D'),fill='#FF950E'),
                    upset_query(set='snp_A', fill='black'),
                    upset_query(set='snp_B', fill='black'),
                    upset_query(set='snp_D', fill='black'),
                    upset_query(set='gene_A', fill='black'),
                    upset_query(set='gene_B', fill='black'),
                    upset_query(set='gene_D', fill='black')
                    
                  )) + xlab("") + theme(axis.text = element_text(size = 18))
transspw

transspwtop <- upset(PWtranssame_upset,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
                     sort_intersections=F,
                     intersections=list(c('snp_A','gene_A'),c('snp_B','gene_B'),c('snp_D','gene_D')),
                     stripes = c(alpha("#FF950E",0.9),alpha("#4B1F6F",0.9),alpha("#579D1C",0.9)),
                     matrix=intersection_matrix(geom = geom_point(size=7),segment = geom_segment(lwd=4)),
                     base_annotations = list(
                       'topsize'=(intersection_size(width=0.3,text = list(size=9),bar_number_threshold = 1) + 
                                    theme_classic() +
                                    theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                          axis.title = element_text(size = 23),
                                          title = element_text(size = 18)) +
                                    coord_cartesian(ylim = c(20000,55000)) +
                                    ylab('') + xlab("") +
                                    ggtitle('PW trans-eQTL same chromosome associations'))
                     ),
                     width_ratio = 0.2,queries = list(
                       upset_query(intersect=c('snp_A','gene_A'),fill='#579D1C'),
                       upset_query(intersect=c('snp_B','gene_B'),fill='#4B1F6F'),
                       upset_query(intersect=c('snp_D','gene_D'),fill='#FF950E'),
                       upset_query(set='snp_A', fill='black'),
                       upset_query(set='snp_B', fill='black'),
                       upset_query(set='snp_D', fill='black'),
                       upset_query(set='gene_A', fill='black'),
                       upset_query(set='gene_B', fill='black'),
                       upset_query(set='gene_D', fill='black')
                       
                     )) + xlab("") + theme(axis.text = element_text(size = 18))
transspwtop

transdpw <- upset(PWtransdiff_upset,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
                  sort_intersections=F,
                  intersections=list(c('snp_A','gene_A'),c('snp_A','gene_B'),c('snp_A','gene_D'),
                                     c('snp_B','gene_A'),c('snp_B','gene_B'),c('snp_B','gene_D'),
                                     c('snp_D','gene_A'),c('snp_D','gene_B'),c('snp_D','gene_D')),
                  stripes = c(alpha("#FF950E",0.9),alpha("#4B1F6F",0.9),alpha("#579D1C",0.9)),
                  matrix=intersection_matrix(geom = geom_point(size=7),segment = geom_segment(lwd=4)),
                  base_annotations = list('size'=(intersection_size(width=0.9,text = list(size=9),bar_number_threshold = 1) + 
                                                    theme_classic() +
                                                    theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                                          axis.title = element_text(size = 23),
                                                          title = element_text(size = 18)) +
                                                    ylab('Number of associations') + xlab("") +
                                                    ylim(0,10000) +
                                                    ggtitle('PW trans-eQTL different chromosome associations'))),
                  width_ratio = 0.2,queries = list(
                    upset_query(intersect=c('snp_A','gene_A'),fill='#579D1C'),
                    upset_query(intersect=c('snp_B','gene_B'),fill='#4B1F6F'),
                    upset_query(intersect=c('snp_D','gene_D'),fill='#FF950E'),
                    upset_query(intersect=c('snp_A','gene_D'),fill='black'),
                    upset_query(intersect=c('snp_B','gene_D'),fill='black'),
                    upset_query(intersect=c('snp_A','gene_B'),fill='black'),
                    upset_query(intersect=c('snp_B','gene_A'),fill='black'),
                    upset_query(intersect=c('snp_D','gene_A'),fill='black'),
                    upset_query(intersect=c('snp_D','gene_B'),fill='black')
                    
                  )) + xlab("") + theme(axis.text = element_text(size = 18))
transdpw

ggsave(plot = transdpw,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transdpw_upset_associations.png',device = 'png',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = transspw,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transspw_upset_associations.png',device = 'png',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = cispw,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/cispw_upset_associations.png',device = 'png',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = transspwtop,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transspwtop_upset_associations.png',device = 'png',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = transdpw,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transdpw_upset_associations.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = transspw,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transspw_upset_associations.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = cispw,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/cispw_upset_associations.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)
ggsave(plot = transspwtop,filename = 'Detection_of_eQTL/new_prep_data/figure_updates/transspwtop_upset_associations.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)






# For Figure 4 barplot of SNPs identified in both populations for each subgenome

snp_postlmmPC <- read.table('Detection_of_eQTL/new_prep_data/postlmm_genotype_PCsinglegene',header = T)
snp_postlmmPW <- read.table('Detection_of_eQTL/new_prep_data/postlmm_genotype_PWsinglegene',header = T)

filteQTLPC <- read.table('Detection_of_eQTL/gatknewway_filteQTLPCC3',header = T)
filteQTLPC <- filteQTLPC[filteQTLPC$id %in% unlist(snp_postlmmPC[,c('id')]),]

filteQTLPW <- read.table('Detection_of_eQTL/gatknewway_filteQTLPW_batch1',header = T)
filteQTLPW <- filteQTLPW[filteQTLPW$id %in% unlist(snp_postlmmPW[,c('id')]),]



pcsnps <- filteQTLPC[,c('chr','gene.id','pos')]
pcsnps$Cross <- 'PxC'
pcsnps$Subgenome <- substr(pcsnps$chr,2,2)

pwsnps <- filteQTLPW[,c('chr','gene.id','pos')]
pwsnps$Cross <- 'PxW'
pwsnps$Subgenome <- substr(pwsnps$chr,2,2)

allsnps <- rbind(pcsnps,pwsnps)
allsnps$Cross <- factor(allsnps$Cross,levels=c('PxC','PxW'))
ggplot(allsnps,aes(y = forcats::fct_rev(Cross),fill=Subgenome)) + geom_bar(position = position_dodge2(reverse = T),color='black') +
  scale_fill_manual(values = c('#579D1C','#4B1F6F','#FF950E')) +
  theme_bw() +
  labs(y='Cross',x='SNP count') +
  xlim(0,5000) +
  theme(axis.title = element_text(size = 38),axis.text = element_text(size = 32),
        legend.title = element_text(size = 32),legend.text = element_text(size = 32)) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxC' & allsnps$Subgenome=='A',])[1]),x=3800,y=2.3,size = 12) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxC' & allsnps$Subgenome=='B',])[1]),x=4400,y=2.0,size = 12) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxC' & allsnps$Subgenome=='D',])[1]),x=1200,y=1.7,size = 12) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxW' & allsnps$Subgenome=='A',])[1]),x=3500,y=1.3,size = 12) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxW' & allsnps$Subgenome=='B',])[1]),x=2400,y=1.0,size = 12) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxW' & allsnps$Subgenome=='D',])[1]),x=1000,y=0.7,size = 12)


ggsave(plot = last_plot(),filename = 'Detection_of_eQTL/new_prep_data/snp_bars.png',device = 'png',width = 12,height = 6,units = 'in',dpi = 500)
ggsave(plot = last_plot(),filename = 'Detection_of_eQTL/new_prep_data/snp_bars.svg',device = 'svg',width = 12,height = 6,units = 'in',dpi = 500)



###### plot the SNP distribution over all chromosomes
chr_len <- read.table('snp_distr_chr_len')
chr_len <- chr_len[1:21,]


pcsnps <- read.table('snp_distr_PC',header = T)
pwsnps <- read.table('snp_distr_PW',header = T)

ggplot() +
  geom_blank(data=chr_len) +
  geom_histogram(data=pcsnps,aes(x=pos,fill = Subgenome,group = chr), position = position_dodge(),binwidth = 8000000) +
  facet_wrap(~chr,ncol = 3,scales = 'free') +
  scale_fill_manual(values = c('#579D1C','#4B1F6F','#FF950E')) +
  theme_bw() +
  labs(y='SNP abundance',x='Position on chromosome')

ggsave('supplementary_figure_SNP_distributionPC.png',device = 'png',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
ggsave('supplementary_figure_SNP_distributionPC.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

ggplot() +
  geom_blank(data=chr_len) +
  geom_histogram(data=pwsnps,aes(x=pos,fill = Subgenome,group = chr), position = position_dodge(),binwidth = 8000000) +
  facet_wrap(~chr,ncol = 3,scales = 'free') +
  scale_fill_manual(values = c('#579D1C','#4B1F6F','#FF950E')) +
  theme_bw() +
  labs(y='SNP abundance',x='Position on chromosome')

ggsave('supplementary_figure_SNP_distributionPW.png',device = 'png',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
ggsave('supplementary_figure_SNP_distributionPW.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
