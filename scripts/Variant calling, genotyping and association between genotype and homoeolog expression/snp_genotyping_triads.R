# This script filters the SNPs and does the haplotyping of homoeologs for the subsequent analysis
library(VariantAnnotation)
library(snpsettest)
library(tidyverse)
library(ggplot2)
library(palettetown)
library(pheatmap)
library(gridExtra)
library(VennDiagram)
library(RColorBrewer)
library(ggtern)
library(randomcoloR)
library(rdist)

# Load in data for genotype and depth per sample
all_vcf <- readVcf("gatk/all3-final.vcf",param=ScanVcfParam(geno = "GT"))
all_vcfGT <- isSNV(all_vcf,singleAltOnly=T)
all_vcfGTdf <- as.data.frame(head(geno(all_vcf))[["GT"]])
all_vcfsnvGT <- subset(all_vcfGTdf,all_vcfGT)
all_vcfdp <- readVcf("gatk/all3-final.vcf",param=ScanVcfParam(geno = "DP"))
all_vcfdp2 <- as.data.frame(head(geno(all_vcfdp))[["DP"]])
all_vcfsnvdp <- subset(all_vcfdp2,all_vcfGT)
colnames(all_vcfsnvGT) <- gsub("-","_",colnames(all_vcfsnvGT))
colnames(all_vcfsnvdp) <- gsub("-","_",colnames(all_vcfsnvdp))
colnames(all_vcfsnvdp) <- paste(colnames(all_vcfsnvdp),"DP",sep="_")
all_vcfsnvdp[is.na(all_vcfsnvdp)] <- 0

rm(all_vcfGTdf)
rm(all_vcfGT)
rm(all_vcfdp2)
rm(all_vcfdp)
# Put them together
###### Free unused memory first!!!
all_vcf <- cbind(all_vcfsnvGT,all_vcfsnvdp)
pcvcfcols <- colnames(all_vcfsnvGT)

all_vcf$split <- row.names(all_vcf)
all_vcf <- separate(data = all_vcf,col = split, into = c("chr","REST"),sep = ":")
all_vcf <- separate(data = all_vcf,col = REST, into = c("pos","SNP"),sep = "_")
all_vcf <- separate(data = all_vcf,col = SNP, into = c("REF","ALT"),sep = "/")
rownames(all_vcf) = NULL
# all_vcf$chr <- gsub("chr","",as.character(all_vcf$chr))
all_vcf$chr = substr(all_vcf$chr,4,5)
all_vcf <- all_vcf[all_vcf$chr != "Un",]

# Remove indels by not allowing more then 1 character in the SNP column
# all_vcf <- all_vcf %>%
#   filter(nchar(REF)<2)


# For snptestset need to create an "snp id column"
all_vcf$id <- 1:nrow(all_vcf)
all_vcf$pos <- as.numeric(all_vcf$pos)

# 
# 
# 
# write out all_vcf preprocessed (CHECKPOINT  to save time when rerunning code)

write.table(all_vcf,file = "/Users/glombik/work/vcf_retry/gatk/recallall_vcf_preprocessed.vcf",quote = F,row.names = F,sep = "\t")
all_vcf <- read.table("/Users/glombik/work/vcf_retry/gatk/recallall_vcf_preprocessed.vcf",header = T)

# And connect it to the annotation file to prepare the base map for viewing
TraeCSgff <- read.table("/Users/glombik/work/vcf_retry/geneIWGSC_v1.1_HC_20170706.gff3sorted")
head(TraeCSgff)
colnames(TraeCSgff) <- c("chr","start","end","gene.id","strand")
TraeCSgff$chr = substr(TraeCSgff$chr,4,5)
TraeCSgff <- TraeCSgff[TraeCSgff$chr != "Un",]

# Now map the SNPs to the annotation
map_all <- map_snp_to_gene(all_vcf,TraeCSgff,extend_start = 0,extend_end = 0)
map_all <- as.data.frame(map_all[["map"]])


map_all$id <- as.numeric(map_all$id)
map_all <- map_all %>%
  dplyr::select(!(gene.adj.start:gene.adj.end)) %>%
  inner_join(all_vcf)

# Check if any NAs - not assigned SNP positions
NAmap_all <- map_all[!complete.cases(map_all),]

# Remove NAs - keep only assigned SNPs
map_all <- map_all[complete.cases(map_all),]
# How many are gene unique?
uniqids <- map_all %>%
  dplyr::count(id)
uniqids_one <- uniqids[uniqids$n == 1,]
uniqids_two <- uniqids[uniqids$n == 2,]
uniqids_thr <- uniqids[uniqids$n == 3,]
uniqids_fr <- uniqids[uniqids$n == 4,]
# where does the SNP belong - take only SNPs that do not overlap between annotated genes
uni_map_all <- map_all %>%
  inner_join(uniqids_one)
uni_map_all$id <- as.numeric(uni_map_all$id)
map_all <- uni_map_all %>%
  dplyr::select(!(n)) %>%
  inner_join(all_vcf)

# # # # # # # # # # # # # # Now, based on the whole genes-with-snp distribution map, extract data for PC cross
# SNPs were already quality filtered during variant calling, do some post-filtering
# Extract only SNPs that are called at least in 10 samples with DP >= 5
sub_all_head <- map_all[,c("id","chr","pos","gene.id","gene.start","gene.end","REF","ALT")]
sub_all_PC <- map_all %>%
  dplyr::select(matches("PC"))
sub_all_PC <- cbind(sub_all_head,sub_all_PC)
qual_sub_all_PC <- sub_all_PC[rowSums(sub_all_PC[9:60] == ".") <=47 & apply(sub_all_PC[,61:112], 1,
                                                                            function(x) length(which(x>=10))) >=5,]

# extract only homozygous SNPs, and filter out uncalled regions
replace_pipe_with_slash <- function(x) {
  gsub("\\|", "/", x)
}
qual_sub_all_PC <- as.data.frame(lapply(qual_sub_all_PC, replace_pipe_with_slash))

reduce_sub_all_PC <- qual_sub_all_PC[with(qual_sub_all_PC, !grepl("0/1", paste(PC_C_,PC_P_))),]
reduce_sub_all_PC <- reduce_sub_all_PC[with(reduce_sub_all_PC, !grepl("\\.", paste(PC_C_,PC_P_))),]
reduce_sub_all_PC <- reduce_sub_all_PC[with(reduce_sub_all_PC, !grepl("1/2", paste(PC_C_,PC_P_))),]

# 268 019 SNPs left
# extract only SNPs that differ between these two parents = can be used to discriminate the homoeolog expression
diffmap_PC <- reduce_sub_all_PC[reduce_sub_all_PC$PC_C_ != reduce_sub_all_PC$PC_P_,]
# 58 845 SNPs left

# # # # # # # # # # # # # # Now, based on the whole genes-with-snp distribution map, extract data for PW cross
# SNPs were already quality filtered during variant calling, do some post-filtering
# Extract only SNPs that are called at least in 10 samples with DP >= 5
sub_all_head <- map_all[,c("id","chr","pos","gene.id","gene.start","gene.end","REF","ALT")]
sub_all_PW <- map_all %>%
  dplyr::select(matches("PW"))
sub_all_PW <- cbind(sub_all_head,sub_all_PW)
qual_sub_all_PW <- sub_all_PW[rowSums(sub_all_PW[9:60] == ".") <=47 & apply(sub_all_PW[,61:112], 1,function(x) length(which(x>=10))) >=5,]


# extract only homozygous SNPs, and filter out uncalled regions
qual_sub_all_PW <- as.data.frame(lapply(qual_sub_all_PW, replace_pipe_with_slash))


reduce_sub_all_PW <- qual_sub_all_PW[with(qual_sub_all_PW, !grepl("0/1", paste(PW_W_,PW_P_))),]
reduce_sub_all_PW <- reduce_sub_all_PW[with(reduce_sub_all_PW, !grepl("\\.", paste(PW_W_,PW_P_))),]
reduce_sub_all_PW <- reduce_sub_all_PW[with(reduce_sub_all_PW, !grepl("1/2", paste(PW_W_,PW_P_))),]

# 259 983 SNPs left
# extract only SNPs that differ between these two parents = can be used to discriminate the homoeolog expression
diffmap_PW <- reduce_sub_all_PW[reduce_sub_all_PW$PW_W_ != reduce_sub_all_PW$PW_P_,]
# 29 262 SNPs left

# Write out diff files
write.table(diffmap_PC,file = "gatk/gatk_diffmap_PC.tab",row.names = F,quote = F,sep = "\t")
write.table(diffmap_PW,file = "gatk/gatk_diffmap_PW.tab",row.names = F,quote = F,sep = "\t")

# Make a barplot of SNPs identified in both populations for each subgenome

pcsnps <- diffmap_PC[,c('chr','gene.id')]
pcsnps$Cross <- 'PxC'
pcsnps$Subgenome <- substr(pcsnps$chr,2,2)

pwsnps <- diffmap_PW[,c('chr','gene.id')]
pwsnps$Cross <- 'PxW'
pwsnps$Subgenome <- substr(pwsnps$chr,2,2)

allsnps <- rbind(pcsnps,pwsnps)
allsnps$Cross <- factor(allsnps$Cross,levels=c('PxC','PxW'))
ggplot(allsnps,aes(y = forcats::fct_rev(Cross),fill=Subgenome)) + geom_bar(position = position_dodge2(reverse = T),color='black') +
  scale_fill_manual(values = c('#579D1C','#4B1F6F','#FF950E')) +
  theme_bw() +
  labs(y='Population',x='SNP count') +
  xlim(0,35000) +
  theme(axis.title = element_text(size = 38),axis.text = element_text(size = 32),
        legend.title = element_text(size = 32),legend.text = element_text(size = 32)) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxC' & allsnps$Subgenome=='A',])[1]),x=30000,y=2.3,size = 12) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxC' & allsnps$Subgenome=='B',])[1]),x=32000,y=2.0,size = 12) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxC' & allsnps$Subgenome=='D',])[1]),x=14000,y=1.7,size = 12) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxW' & allsnps$Subgenome=='A',])[1]),x=21000,y=1.3,size = 12) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxW' & allsnps$Subgenome=='B',])[1]),x=15000,y=1.0,size = 12) +
  annotate("text",label = paste0(dim(allsnps[allsnps$Cross=='PxW' & allsnps$Subgenome=='D',])[1]),x=10000,y=0.7,size = 12)


ggsave(plot = last_plot(),filename = 'snp_bars.png',device = 'png',width = 12,height = 6,units = 'in',dpi = 500)
ggsave(plot = last_plot(),filename = 'snp_bars.svg',device = 'svg',width = 12,height = 6,units = 'in',dpi = 500)

table(allsnps$Cross,allsnps$Subgenome)



# This script is gonna compare the SNPs with the expression bias


# load in the SNP parental-differential data
diffmap_PC <- read.table("gatk/gatk_diffmap_PC.tab",header = T)
diffmap_PW <- read.table("gatk/gatk_diffmap_PW.tab",header = T)

#### Check how similar they are since the numbers are pretty much the same
freebdiffmapPC <- read.table("diffmap_PC.tab",header=T)
freebdiffmapPW <- read.table("diffmap_PW.tab",header=T)
freebdiffmapPC <- freebdiffmapPC[,2:7]
freebdiffmapPW <- freebdiffmapPW[,2:7]
freebdiffmapPC$pos <- as.character(freebdiffmapPC$pos)
freebdiffmapPW$pos  <- as.character(freebdiffmapPW$pos)

diffmap_PC$gene.start <- as.integer(diffmap_PC$gene.start)
diffmap_PC$gene.end <- as.integer(diffmap_PC$gene.end)

PCdiffcomp <- diffmap_PC %>%
  inner_join(freebdiffmapPC)
# 50 138

diffmap_PW$gene.start <- as.integer(diffmap_PW$gene.start)
diffmap_PW$gene.end <- as.integer(diffmap_PW$gene.end)

PWdiffcomp <- diffmap_PW %>%
  inner_join(freebdiffmapPW)
# 24 344
##########

# FIRST WE WILL FILTER THE SNP DATA, TRANSFER IT INTO LONG FORMAT AND CREATE A COLUMN FOR EACH F5 SAMPLE THAT SAYS IF
# THE GENOTYPE CALLED IS: PARENT1 (P1,Paragon) / PARENT2 (P2,Charger/Watkins) / NOT-CALLED (NC) / HETEROZYGOUS (HET)

# get rid of depth columns
diffmap_PCnoDP <- diffmap_PC %>%
  select(!matches("DP"))
# take only F5 samples (no parents)
f6diffmap_PCnoDP <- diffmap_PCnoDP %>%
  select(!matches(c("PC_C","PC_P")))
# take only parents here
pardiffmap_PCnoDP <- diffmap_PCnoDP %>%
  select(matches(c("id","chr","pos","gene.id","gene.start","gene.end","REF","ALT","PC_C","PC_P")))
# make f5 into long format
lngdiffmapPC <- f6diffmap_PCnoDP %>% pivot_longer(!id:ALT,names_to = "sample")
colnames(lngdiffmapPC)[10] <- c("genotype")
# now join long format f5 back to wide format parents
combined_diff_PC <- lngdiffmapPC %>%
  inner_join(pardiffmap_PCnoDP)
combined_diff_PC$GT_match <- "NC"
combined_diff_PC$GT_match[combined_diff_PC$genotype == combined_diff_PC$PC_P_] <- "P1"
combined_diff_PC$GT_match[combined_diff_PC$genotype == combined_diff_PC$PC_C_] <- "P2"
combined_diff_PC$GT_match[combined_diff_PC$genotype == "0/1"] <- "HET"
combined_diff_PC$GT_match[combined_diff_PC$genotype == "1/2"] <- "HET"



diffmap_PWnoDP <- diffmap_PW %>%
  select(!matches("DP"))
# take only F5 samples (no parents)
f6diffmap_PWnoDP <- diffmap_PWnoDP %>%
  select(!matches(c("PW_W","PW_P")))
# take only parents here
pardiffmap_PWnoDP <- diffmap_PWnoDP %>%
  select(matches(c("id","chr","pos","gene.id","gene.start","gene.end","REF","ALT","PW_W","PW_P")))
# make into long format
lngdiffmapPW <- f6diffmap_PWnoDP %>% pivot_longer(!id:ALT,names_to = "sample")
colnames(lngdiffmapPW)[10] <- c("genotype")
# now join long format f5 back to wide format parents
combined_diff_PW <- lngdiffmapPW %>%
  inner_join(pardiffmap_PWnoDP)
combined_diff_PW$GT_match <- "NC"
combined_diff_PW$GT_match[combined_diff_PW$genotype == combined_diff_PW$PW_P_] <- "P1"
combined_diff_PW$GT_match[combined_diff_PW$genotype == combined_diff_PW$PW_W_] <- "P2"
combined_diff_PW$GT_match[combined_diff_PW$genotype == "0/1"] <- "HET"
combined_diff_PW$GT_match[combined_diff_PW$genotype == "1/2"] <- "HET"

# Haplotype the genes according to SNPs?

# Haplotype test ----------------------------------------------------------

# Check how few samples look before running the haplotyping...PC_10, PW_8

pc_10 <- combined_diff_PC[combined_diff_PC$sample == "PC_10",]

# Need to count gene.id+GT_match to see how many genes show nice 1 GT pattern and which are problematic
PCuniq_geneGTs <- pc_10 %>%
  dplyr::count(gene.id,GT_match)
PCuniq_gene <- pc_10 %>%
  dplyr::count(gene.id)
colnames(PCuniq_geneGTs)[3] <- "nGT"
colnames(PCuniq_gene)[2] <- "ngene"
PCuniqmerge <- PCuniq_geneGTs %>%
  inner_join(PCuniq_gene)
PCuniqwide <- PCuniqmerge %>%
  pivot_wider(names_from = GT_match,values_from = nGT)

PCuniqwide[is.na(PCuniqwide)] <- 0
PCuniqwide$HET.ratio <- PCuniqwide$HET/PCuniqwide$ngene
PCuniqwide$P2.ratio <- PCuniqwide$P2/PCuniqwide$ngene
PCuniqwide$P1.ratio <- PCuniqwide$P1/PCuniqwide$ngene
PCuniqwide$NC.ratio <- PCuniqwide$NC/PCuniqwide$ngene
P2full <- PCuniqwide[PCuniqwide$P2.ratio==1,]
P1full <- PCuniqwide[PCuniqwide$P1.ratio==1,]
HETfull <- PCuniqwide[PCuniqwide$HET.ratio==1,]
NCfull <- PCuniqwide[PCuniqwide$NC.ratio==1,]
# P2 = 1...5081, P1 = 1...4135, HET = 847, NC = 453, total = 12325

P2fullorplusNC <- PCuniqwide[PCuniqwide$P2.ratio==1 | 
                               (PCuniqwide$P2.ratio + PCuniqwide$NC.ratio == 1 & 
                                  PCuniqwide$NC.ratio > 0 & PCuniqwide$P2.ratio > 0),]
P1fullorplusNC <- PCuniqwide[PCuniqwide$P1.ratio==1 |
                               (PCuniqwide$P1.ratio + PCuniqwide$NC.ratio == 1 & 
                                  PCuniqwide$NC.ratio > 0 & PCuniqwide$P1.ratio > 0),]
HETfullorplusNC <- PCuniqwide[PCuniqwide$HET.ratio==1 |
                                (PCuniqwide$HET.ratio + PCuniqwide$NC.ratio == 1 & 
                                   PCuniqwide$NC.ratio > 0 & PCuniqwide$HET.ratio > 0),]

# Plot the SNP density per gene
g1 <- ggplot(PCuniqwide) + geom_bar(aes(ngene),fill="pink") +
  xlab("SNPs per gene") +
  ylab("Count") +
  geom_rect(aes(xmin=0,xmax=20,ymin=0,ymax=4700),color="black", alpha=0) +
  theme_classic()

g2 <- ggplot(PCuniqwide) + geom_bar(aes(ngene),fill="pink") + xlim(0,20) +
  xlab("SNPs per gene") +
  ylab("Count") +
  theme_classic()

g12 <- g1 + annotation_custom(ggplotGrob(g2),xmin = 25, xmax = 75,ymin = 1000,ymax = 4000) +
  geom_rect(aes(xmin=25,xmax=75,ymin=1000,ymax=4000),color="black", linetype='dashed', alpha=0) +
  geom_path(aes(x,y,group=grp),data=data.frame(x=c(20,25,20,25),y=c(0,1000,4700,4000),grp=c(1,1,2,2)),
            linetype="dashed") +
  ggtitle("PxC")

# ggsave(filename = "SNPpergene_dens_PC.png",device = "png",dpi = 400,plot = g12,width = 17,height = 10,units = "in")


# Now the PW sample

PW_8 <- combined_diff_PW[combined_diff_PW$sample == "PW_8_",]

# Need to count gene.id+GT_match to see how many genes show nice 1 GT pattern and which are problematic
PWuniq_geneGTs <- PW_8 %>%
  dplyr::count(gene.id,GT_match)
PWuniq_gene <- PW_8 %>%
  dplyr::count(gene.id)
colnames(PWuniq_geneGTs)[3] <- "nGT"
colnames(PWuniq_gene)[2] <- "ngene"
PWuniqmerge <- PWuniq_geneGTs %>%
  inner_join(PWuniq_gene)
PWuniqwide <- PWuniqmerge %>%
  pivot_wider(names_from = GT_match,values_from = nGT)

PWuniqwide[is.na(PWuniqwide)] <- 0
PWuniqwide$HET.ratio <- PWuniqwide$HET/PWuniqwide$ngene
PWuniqwide$P2.ratio <- PWuniqwide$P2/PWuniqwide$ngene
PWuniqwide$P1.ratio <- PWuniqwide$P1/PWuniqwide$ngene
PWuniqwide$NC.ratio <- PWuniqwide$NC/PWuniqwide$ngene
P2full <- PWuniqwide[PWuniqwide$P2.ratio==1,]
P1full <- PWuniqwide[PWuniqwide$P1.ratio==1,]
HETfull <- PWuniqwide[PWuniqwide$HET.ratio==1,]
NCfull <- PWuniqwide[PWuniqwide$NC.ratio==1,]

P2fullorplusNC <- PWuniqwide[PWuniqwide$P2.ratio==1 | 
                               (PWuniqwide$P2.ratio + PWuniqwide$NC.ratio == 1 & 
                                  PWuniqwide$NC.ratio > 0 & PWuniqwide$P2.ratio > 0),]
P1fullorplusNC <- PWuniqwide[PWuniqwide$P1.ratio==1 |
                               (PWuniqwide$P1.ratio + PWuniqwide$NC.ratio == 1 & 
                                  PWuniqwide$NC.ratio > 0 & PWuniqwide$P1.ratio > 0),]
HETfullorplusNC <- PWuniqwide[PWuniqwide$HET.ratio==1 |
                                (PWuniqwide$HET.ratio + PWuniqwide$NC.ratio == 1 & 
                                   PWuniqwide$NC.ratio > 0 & PWuniqwide$HET.ratio > 0),]




# Haplotyping -------------------------------------------------------------
# PC Haplotyping
PCuniq_geneGTs <- combined_diff_PC %>%
  dplyr::count(gene.id,GT_match,sample)
PCuniq_gene <- combined_diff_PC %>%
  dplyr::count(gene.id,sample)
colnames(PCuniq_geneGTs)[4] <- "nGT"
colnames(PCuniq_gene)[3] <- "ngene"
PCuniqmerge <- PCuniq_geneGTs %>%
  inner_join(PCuniq_gene)
PCuniqwide <- PCuniqmerge %>%
  pivot_wider(names_from = GT_match,values_from = nGT)

PCuniqwide[is.na(PCuniqwide)] <- 0
PCuniqwide$HET.ratio <- PCuniqwide$HET/PCuniqwide$ngene
PCuniqwide$P2.ratio <- PCuniqwide$P2/PCuniqwide$ngene
PCuniqwide$P1.ratio <- PCuniqwide$P1/PCuniqwide$ngene
PCuniqwide$NC.ratio <- PCuniqwide$NC/PCuniqwide$ngene


PCuniqwide$GTfilter <- "ND" # Not Distinguishable
PCuniqwide$GTfilter[PCuniqwide$P2.ratio== 1 | (PCuniqwide$P2.ratio + PCuniqwide$NC.ratio == 1 & PCuniqwide$NC.ratio > 0 & PCuniqwide$P2.ratio > 0)] <- "P2"
PCuniqwide$GTfilter[PCuniqwide$P1.ratio== 1 | (PCuniqwide$P1.ratio + PCuniqwide$NC.ratio == 1 & PCuniqwide$NC.ratio > 0 & PCuniqwide$P1.ratio > 0)] <- "P1"
PCuniqwide$GTfilter[PCuniqwide$HET.ratio== 1 | (PCuniqwide$HET.ratio + PCuniqwide$NC.ratio == 1 & PCuniqwide$NC.ratio > 0 & PCuniqwide$HET.ratio > 0)] <- "HET"

PCuniqwide %>%
  count(sample,GTfilter) %>%
  pivot_wider(names_from = GTfilter,values_from = n) %>%
  summarise(mean_HET = mean(HET),
            mean_P1 = mean(P1),
            mean_P2 = mean(P2),
            mean_ND = mean(ND))

PC_GTfilter_plot <- ggplot(PCuniqwide) + geom_bar(aes(x=sample,fill=GTfilter),alpha=0.7) +
  annotate(geom="rect",xmin=29.5,xmax=30.5,ymin=0,ymax=12325,color="red", alpha=0) +
  annotate(geom="rect",xmin=48.5,xmax=49.5,ymin=0,ymax=12325,color="red", alpha=0) +
  geom_hline(yintercept = 3000) +
  geom_hline(yintercept = 8000) +
  geom_hline(yintercept = 11000) +
  geom_hline(yintercept = 12000) +
  annotate("text",label = "mean P1 = 4912",x=5,y=2800,size = 4) +
  annotate("text",label = "mean P2 = 4985",x=5,y=7800,size = 4) +
  annotate("text",label = "mean ND = 1367",x=5,y=10800,size = 4) +
  annotate("text",label = "mean HET = 1061",x=5,y=11800,size = 4) +
  ylab("count") +
  scale_fill_poke(pokemon = "Mewtwo") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

ggsave(filename = "gatk/PC_GTfilter_plot.png",device = "png",dpi = 400,plot = PC_GTfilter_plot,width = 10,height = 10,units = "in")

# PW Haplotyping
PWuniq_geneGTs <- combined_diff_PW %>%
  dplyr::count(gene.id,GT_match,sample)
PWuniq_gene <- combined_diff_PW %>%
  dplyr::count(gene.id,sample)
colnames(PWuniq_geneGTs)[4] <- "nGT"
colnames(PWuniq_gene)[3] <- "ngene"
PWuniqmerge <- PWuniq_geneGTs %>%
  inner_join(PWuniq_gene)
PWuniqwide <- PWuniqmerge %>%
  pivot_wider(names_from = GT_match,values_from = nGT)

PWuniqwide[is.na(PWuniqwide)] <- 0
PWuniqwide$HET.ratio <- PWuniqwide$HET/PWuniqwide$ngene
PWuniqwide$P2.ratio <- PWuniqwide$P2/PWuniqwide$ngene
PWuniqwide$P1.ratio <- PWuniqwide$P1/PWuniqwide$ngene
PWuniqwide$NC.ratio <- PWuniqwide$NC/PWuniqwide$ngene

PWuniqwide$GTfilter <- "ND" # Not Distinguishable
PWuniqwide$GTfilter[PWuniqwide$P2.ratio== 1 | (PWuniqwide$P2.ratio + PWuniqwide$NC.ratio == 1 & PWuniqwide$NC.ratio > 0 & PWuniqwide$P2.ratio > 0)] <- "P2"
PWuniqwide$GTfilter[PWuniqwide$P1.ratio== 1 | (PWuniqwide$P1.ratio + PWuniqwide$NC.ratio == 1 & PWuniqwide$NC.ratio > 0 & PWuniqwide$P1.ratio > 0)] <- "P1"
PWuniqwide$GTfilter[PWuniqwide$HET.ratio== 1 | (PWuniqwide$HET.ratio + PWuniqwide$NC.ratio == 1 & PWuniqwide$NC.ratio > 0 & PWuniqwide$HET.ratio > 0)] <- "HET"

PWuniqwide %>%
  count(sample,GTfilter) %>%
  pivot_wider(names_from = GTfilter,values_from = n) %>%
  summarise(mean_HET = mean(HET),
            mean_P1 = mean(P1),
            mean_P2 = mean(P2),
            mean_ND = mean(ND))

PW_GTfilter_plot <-ggplot(PWuniqwide) + geom_bar(aes(x=sample,fill=GTfilter),alpha=0.7) +
  annotate(geom="rect",xmin=4.5,xmax=5.5,ymin=0,ymax=6447,color="red", alpha=0) +
  annotate(geom="rect",xmin=15.5,xmax=16.5,ymin=0,ymax=6447,color="red", alpha=0) +
  annotate(geom="rect",xmin=28.5,xmax=29.5,ymin=0,ymax=6447,color="red", alpha=0) +
  annotate(geom="rect",xmin=38.5,xmax=39.5,ymin=0,ymax=6447,color="red", alpha=0) +
  geom_hline(yintercept = 1400) +
  geom_hline(yintercept = 4000) +
  geom_hline(yintercept = 5800) +
  geom_hline(yintercept = 6300) +
  annotate("text",label = "mean P1 = 2639",x=5,y=1200,size = 4) +
  annotate("text",label = "mean P2 = 2495",x=5,y=3800,size = 4) +
  annotate("text",label = "mean ND = 762",x=5,y=5600,size = 4) +
  annotate("text",label = "mean HET = 551",x=5,y=6100,size = 4) +
  ylab("count") +
  scale_fill_poke(pokemon = "Mewtwo") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1))

ggsave(filename = "gatk/PW_GTfilter_plot.png",device = "png",dpi = 400,plot = PW_GTfilter_plot,width = 10,height = 10,units = "in")



# Merge with homologies and triads ----------------------------------------

# load in the homology data
homologies <- read.csv(file="/Users/glombik/work/obj1_reanalysis/homoeologs_1_1_1_synt_and_non_synt.csv")

lnghomologies <- homologies %>%
  pivot_longer(A:D,names_to = "genome")
colnames(lnghomologies)[8] <- c("gene.id")

PCdiff_homol <- PCuniqwide %>%
  inner_join(lnghomologies) %>%
  dplyr::select(!(ngene:NC.ratio))

PCdiff_wide <- PCdiff_homol %>%
  pivot_wider(names_from = genome,values_from = c(gene.id,GTfilter),) %>%
  left_join(homologies) %>%
  dplyr::select(!(gene.id_A:gene.id_D))

PCdiff_wide$chrs[is.na(PCdiff_wide$chrs)] <- as.character("Un")

PCdiff_wide$sample <- gsub("_$","",PCdiff_wide$sample)
colnames(PCdiff_wide)[1] <- "genotype"

write.table(PCdiff_wide,"/Users/glombik/work/vcf_retry/gatk/PCdiff_wide",row.names = F)


PWdiff_homol <- PWuniqwide %>%
  inner_join(lnghomologies) %>%
  dplyr::select(!(ngene:NC.ratio))

PWdiff_wide <- PWdiff_homol %>%
  pivot_wider(names_from = genome,values_from = c(gene.id,GTfilter),) %>%
  left_join(homologies) %>%
  dplyr::select(!(gene.id_A:gene.id_D))

PWdiff_wide$chrs[is.na(PWdiff_wide$chrs)] <- as.character("Un")

PWdiff_wide$sample <- gsub("_$","",PWdiff_wide$sample)
colnames(PWdiff_wide)[1] <- "genotype"

write.table(PWdiff_wide,"/Users/glombik/work/vcf_retry/gatk/PWdiff_wide",row.names = F)




