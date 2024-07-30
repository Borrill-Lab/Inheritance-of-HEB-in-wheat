# Arunkumar Ramesh
# This script calculates CV in triads in F5 lines and computes initial bias distance calculations
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
library(tidyverse)
library(dplyr)
library(pracma)
library(ggplot2)
library(scran)
library(LSD)
library(readr)
library(lme4)
library(ggpubr)
library(car)
library(multcomp)
library(MASS)
library(pheatmap)
library(cowplot)
library(ggvenn)
library(equivalence)
library(BayesFactor)
# library(DEGreport)

TPM <- read.csv(file="bias_category_all_samples_inc_orig_expr.csv")
TPM$sample <- gsub("PC.","PC_",TPM$sample)
TPM$sample <- gsub("PW.","PW_",TPM$sample)
TPM$genotype <- gsub("\\..*","",TPM$sample)

## this part calculates the CV2 and distance to median (DM) for each gene across the three replicates for each genotype, loop takes a long while to finish
TPM_melt <- TPM %>% pivot_longer(cols =c(1:3), names_to = "subgenome", values_to = "TPM")                           
TPM_melt$gene <- paste(TPM_melt$group_id,TPM_melt$subgenome,sep="_")
genotypes <- unique(TPM$genotype)
TPM_g_all <- ""
for (g in 1:length(genotypes)){
  TPM_g <- TPM_melt[TPM_melt$genotype %in% genotypes[g],]
  TPM_g <- TPM_g %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(mean_tpm = mean(TPM), CV2 = ((sd(TPM))^2)/((mean(TPM))^2) )
  TPM_g$DM <- DM(TPM_g$mean_tpm,TPM_g$CV2)
  TPM_g$genotype <- genotypes[g]
  TPM_g_all <- rbind(TPM_g_all,TPM_g)
}

TPM_g_all <- TPM_g_all[-c(1),]
TPM_g_all[2:4] <- apply(TPM_g_all[2:4],2,as.numeric)
write.csv(TPM_g_all,file="TPM_g_all.csv",row.names = F, quote = F)

TPM_g <- read.csv(file="TPM_g_all.csv")
TPM_g <- TPM_g[TPM_g$mean_tpm > 0.5,] ## only analyse genes whose mean TPM across the 3 replicates per genotype is > 0.5 
# here plot distance to median (estimate of expression noise) by genotype and traid
TPM_g2 <- TPM_g %>%
  dplyr::group_by(genotype) %>%
  dplyr::summarise(mean_dm = median(DM,na.rm=T))
TPM_g2 <- TPM_g2[order(TPM_g2$mean_dm),]
ggplot(data=TPM_g,aes(x=factor(genotype,levels=TPM_g2$genotype),y=DM))+
  geom_boxplot()

# for triads, only plotting a random sample of 500
genes <- unique(TPM_g$gene)
TPM_g2 <- TPM_g %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(mean_dm = median(DM,na.rm=T))
TPM_g2 <- TPM_g2[order(TPM_g2$mean_dm),]
ggplot(data=TPM_g[TPM_g$gene %in% sample(genes,500),],aes(x=factor(gene,levels=TPM_g2$gene),y=DM))+
  geom_boxplot()

## this goal of this section is to calculate the area of the triangle (where the three points represent the 3 replicates per genotype) in 3D space (where x,y and z axes represent the A,B and D subgenome expression)
TPM <- TPM[order(TPM$genotype,TPM$group_id),]
a <- TPM %>%
  dplyr::group_by(genotype,group_id) %>%
  dplyr::mutate(A_new = A-dplyr::first(A))
a <- a %>%
  dplyr::group_by(genotype,group_id) %>%
  dplyr::mutate(B_new = B-dplyr::first(B))
a <- a %>%
  dplyr::group_by(genotype,group_id) %>%
  dplyr::mutate(D_new = D-dplyr::first(D))
## only keep triads that occur in all 3 replicates per genotype

a <- a %>%
  dplyr::group_by(genotype,group_id) %>%
  dplyr::filter(n()>2)

a1 <- a[seq(2, nrow(a), 3),((ncol(a))-2):ncol(a)]
a2 <- a[seq(3, nrow(a), 3),((ncol(a))-2):ncol(a)]
alist <- ""
for(l in 1:nrow(a1)){
  print(l)
  alist[l] <- (sqrt(sum(pracma::cross(as.matrix(a1[l,]),as.matrix(a2[l,]))^2)))/2
}

n <- a %>%
  dplyr::group_by(genotype,group_id) %>% 
  dplyr::summarise(n())
n$area <- as.numeric(alist)
n$id <- paste(n$group_id,n$genotype,sep="_")

n2 <- n[!is.na(n$area),]
n2 <- n2[!is.infinite(n2$area),]

write.csv(n2,file="n2.csv",row.names = F)

### this section below explores area of noise by genotype and by triads

n2 <- read_csv("n2.csv")

quantile(-log(n2$area),na.rm=T)
quantile(-log(n2$area),na.rm=T,probs=c(0.01,0.025,0.05))

ggplot(data=n2,aes(x=-log(area)))+
  geom_histogram(bins=100)+
  geom_vline(aes(xintercept = -log(0.01)))

ggplot(data=n2,aes(x=genotype,y=-log(area)))+
  geom_boxplot()+
  geom_hline(aes(yintercept = -log(0.01)))

triadids <- unique(n2$group_id)
n3 <- n2 %>%
  dplyr::group_by(group_id) %>%
  dplyr::summarise(mean_area = median(area,na.rm=T))
n3 <- n3[order(n3$mean_area),]
ggplot(data=n2[n2$group_id %in% sample(triadids,500),],aes(x=factor(group_id,levels=n3$group_id),y=-log(area)))+
  geom_boxplot()

## merging noise area with TPM counts
n <- n2[-c(1:3)]
TPM$id <- paste(TPM$group_id,TPM$genotype,sep="_")
TPM <- left_join(TPM,n,by="id")
TPM <- TPM[order(-TPM$area),]
length(unique(TPM$group_id)) # 16844 triads ---- MAREK UPDATE 16746
TPM <- TPM[-c(7:8)]
#rm(a)

## number of triads by sample
numtraidsbygt <- TPM %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(num_triads = n())
ggplot(numtraidsbygt,aes(x=sample,y=num_triads))+
  geom_col()
numtraidsbygt[numtraidsbygt$num_triads < 12000,]

## contains library read depth estimates from genewiz
seq_stat <- read.csv("40-550500306 Table 2.1 Sample sequencing statistics_first_delivery.csv")
colnames(seq_stat)[2] <- "Genewiz.ID"
seq_stat <- seq_stat[order(seq_stat$Genewiz.ID),]
### note that for 18 samples, a separate library prep and sequencing was done. for the purposes of samples read depth, values are summed but library details differ
seq_stat2 <- read.csv("40-550500306 Table 2.1 Sample sequencing statistics_second_delivery.csv")
colnames(seq_stat2)[2] <- "Genewiz.ID"
seq_stat2 <- seq_stat2[order(seq_stat2$Genewiz.ID),]
seq_stat2$Genewiz.ID <- gsub("ARR02_PC-49-1-2","ARR01_PC-49-1",seq_stat2$Genewiz.ID)
seq_stat <- rbind(seq_stat,seq_stat2[!seq_stat2$Genewiz.ID %in% seq_stat$Genewiz.ID,]) ## last two lines actually for third and fourth delivery now added to first seqstat df
seq_stat2 <- seq_stat2[nrow(seq_stat2),] ## last two lines actually for third and fourth delivery, removed now
seq_stat2 <- rbind(seq_stat2,seq_stat[seq_stat$Genewiz.ID %in% "ARR01_PW-30-2-2",])
seq_stat2$Genewiz.ID <- gsub("ARR01_PW-30-2-2","ARR01_PW-30-2",seq_stat2$Genewiz.ID)
seq_stat <- seq_stat[!seq_stat$Genewiz.ID %in% "ARR01_PW-30-2-2",]
seq_stat[seq_stat$Genewiz.ID %in% seq_stat2$Genewiz.ID,]$Yield..Mbases. <- seq_stat[seq_stat$Genewiz.ID %in% seq_stat2$Genewiz.ID,]$Reads + seq_stat2$Yield..Mbases.
seq_stat[seq_stat$Genewiz.ID %in% seq_stat2$Genewiz.ID,]$Yield..Mbases. <- seq_stat[seq_stat$Genewiz.ID %in% seq_stat2$Genewiz.ID,]$Reads + seq_stat2$Yield..Mbases.
colnames(seq_stat)[2] <- "sample"
seq_stat$sample <- gsub("ARR01_PC-","PC_",seq_stat$sample)
seq_stat$sample <- gsub("ARR01_PW-","PW_",seq_stat$sample)
seq_stat$sample <- gsub("-",".",seq_stat$sample)
seq_stat <- cbind(seq_stat[2],seq_stat[4])

numtraidsbygt <- dplyr::inner_join(numtraidsbygt,seq_stat,by="sample")

## number of triads by number of reads
ggscatter(numtraidsbygt,x="Reads",y="num_triads")+
  xlab("Number of reads")+
  ylab("Number of triads")

## CV (measure of expression bias) and triad tpm
TPM$cv <- apply(TPM[1:3],1,sd)/apply(TPM[1:3],1,mean)
TPM$triad_tpm <- TPM$A_tpm+TPM$B_tpm+TPM$D_tpm


write.table(TPM,file = "TPMbeforeHEB_line193",row.names = F,quote = F,sep = "\t")
# ###### CHECKPOINT - intermediate file for later scripts to use


TPM <- read.table("/Users/glombik/work/obj1_reanalysis/TPMbeforeHEB_line193",header = T)
# can load the file above
## summarize tpm per genotype
HEB <- TPM %>%
  dplyr::group_by(genotype,group_id,id,area) %>% 
  dplyr::summarise(A_tpm=mean(A_tpm),B_tpm=mean(B_tpm),D_tpm=mean(D_tpm))

HEB$A <- HEB$A_tpm/(HEB$A_tpm+HEB$B_tpm+HEB$D_tpm)
HEB$B <- HEB$B_tpm/(HEB$A_tpm+HEB$B_tpm+HEB$D_tpm)
HEB$D <- HEB$D_tpm/(HEB$A_tpm+HEB$B_tpm+HEB$D_tpm)
HEB$group_id <- as.character(HEB$group_id)

## since summarizing by genotype, only keep triads and genotypes with low between replicate variaiton
HEB <- HEB[HEB$area < 0.01,]

## CV (measure of expression bias) and triad tpm

HEB$cv <- apply(HEB[5:7],1,sd)/apply(HEB[5:7],1,mean)
HEB$triad_tpm <- HEB$A_tpm+HEB$B_tpm+HEB$D_tpm

## slight negative association (r=-0.29) between expression bias and expression level
heatscatter(x=log2(HEB$triad_tpm),y=HEB$cv,main="",cor=T)

## first to PC specific analyses
## this part tests for triads with significant evidence of expression bias

## obtain gene expression data, note using all replicates, not just mean per genotype for stats
lmerdata <- TPM[grepl("PC_",TPM$genotype),]
lmerdata$sample <- gsub("\\.","-", lmerdata$sample)
lmerdata$sample <- gsub("-1-2$","-1",lmerdata$sample)
lmerdata$sample <- gsub("-2-2$","-2",lmerdata$sample)
lmerdata$sample <- gsub("-3-2$","-3",lmerdata$sample)
lmerdata$sample <- gsub("-4-2$","-4",lmerdata$sample)
lmerdata$replicate <- gsub(".*-","", lmerdata$sample )
lmerdata$group_id <- as.character(lmerdata$group_id)
rownames(lmerdata) <- paste(lmerdata$group_id,lmerdata$sample,sep=",")

## do a test for cv with genotype as fixed effect and replicate as random effect. 
## only for triads with at least 50 samples assayed
gentoypecvpvals_pc <- lmerdata[!lmerdata$group_id %in% "16241",] %>%
  dplyr::group_by(group_id) %>%
  dplyr::filter(n()>50) %>%
  dplyr::summarise(pval = Anova(lmer(formula=cv~genotype+(1|replicate)))[,3])
## adjust pvalues using fdr to correct for multiple testing
gentoypecvpvals_pc$adjustP <- p.adjust(gentoypecvpvals_pc$pval,method = "BH")
gentoypecvpvals_pc$gt_sig <- T
gentoypecvpvals_pc[gentoypecvpvals_pc$adjustP > 0.05,]$gt_sig <- F
gentoypecvpvals_pc2 <- gentoypecvpvals_pc[-c(2:3)]
table(gentoypecvpvals_pc2$gt_sig)

## note that I also tried to using a negative binomial model (via limma-voom) to test for interaction between genotype and subgenome on TPM (implemented in de.R). the model did not work well when a) CV did not change but TPM did in just one subgenome and b) when variance in genotype much larger than variance in subgenome

# example triads not significant in lmer(cv ~ genotype) but significant in glm.nb(TPM ~ genotype * subgenome)
#"2279","15638","4305", "2051","17304"
pheatmap(t(lmerdata[lmerdata$group_id %in% "15638",4:6]),cluster_cols = F) ## relative expression
pheatmap(t(lmerdata[lmerdata$group_id %in% "15638",1:3]),cluster_cols = F) ## total expression

# example triads  significant in lmer(cv ~ genotype) but not significant in glm.nb(TPM ~ genotype * subgenome)
#"7972","13581","13581","8143","3150"
pheatmap(t(lmerdata[lmerdata$group_id %in% "13581",4:6]),cluster_cols = F)
pheatmap(t(lmerdata[lmerdata$group_id %in% "13581",1:3]),cluster_cols = F)

# example triads not significant in both lmer(cv ~ genotype) and glm.nb(TPM ~ genotype * subgenome)
#"10142", "5625", "2639", "14751", "680"
pheatmap(t(lmerdata[lmerdata$group_id %in% "5625",4:6]),cluster_cols = F)
pheatmap(t(lmerdata[lmerdata$group_id %in% "5625",1:3]),cluster_cols = F)

# what is a good pvalue threshold? Visually, 0.0001 seems like a good threshold but it is subjective
gentoypecvpvals_pc <- gentoypecvpvals_pc[order(gentoypecvpvals_pc$adjustP),]
## p = 0.05
tail(gentoypecvpvals_pc[gentoypecvpvals_pc$adjustP < 0.05,])
pheatmap(t(lmerdata[lmerdata$group_id %in% "6516",4:6]),cluster_cols = F)
## p = 0.01
tail(gentoypecvpvals_pc[gentoypecvpvals_pc$adjustP < 0.01,])
pheatmap(t(lmerdata[lmerdata$group_id %in% "2031",4:6]),cluster_cols = F)
## p = 0.001
tail(gentoypecvpvals_pc[gentoypecvpvals_pc$adjustP < 0.001,])
pheatmap(t(lmerdata[lmerdata$group_id %in% "9913",4:6]),cluster_cols = F)
## p = 0.0001
tail(gentoypecvpvals_pc[gentoypecvpvals_pc$adjustP < 0.0001,])
pheatmap(t(lmerdata[lmerdata$group_id %in% "11618",4:6]),cluster_cols = F)
## p = 1
# tail(gentoypecvpvals_pc)
# pheatmap(t(lmerdata[lmerdata$group_id %in% "11618",4:6]),cluster_cols = F)

## only analyze triads with significant evidence of expression bias variation between genotypes
HEB_PC <- HEB[grepl("PC_",HEB$genotype),]
HEB_PC <- HEB_PC[HEB_PC$group_id %in% gentoypecvpvals_pc[gentoypecvpvals_pc$adjustP < 0.0001,]$group_id,]

## parental triad expression estimates
PC_P <- HEB_PC %>%
  dplyr::filter(genotype == "PC_P")
PC_C <- HEB_PC %>%
  dplyr::filter(genotype == "PC_C")

PC_P_C <- inner_join(PC_P,PC_C,by="group_id")

PC_P_cv_tpm <- cbind(PC_P[c(2,11,12)],"To P")
colnames(PC_P_cv_tpm)[2:4] <- c("cv_parent","tpm_parent", "which_parent")
PC_C_cv_tpm <- cbind(PC_C[c(2,11,12)],"To C")
colnames(PC_C_cv_tpm)[2:4] <- c("cv_parent","tpm_parent", "which_parent")

## calculates distance matrix between parents
## this distance matrix is estimated by comparing first three columns to last three columns for each row separately
PC_all <- cbind(PC_P_C[8:10],PC_P_C[19:21])
PC_P_C$distance <- apply(PC_all, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
PC_P_C_sum <- cbind(PC_P_C[2],PC_P_C[ncol(PC_P_C)])

pdf("PC parents compare.pdf",height=3,width=9)
par(mfrow = c(1, 3))

## compares HEB between parents (r=0.58)
heatscatter(x=PC_P_C$cv.x,y=PC_P_C$cv.y,cor=T,main="",xlab="expression bias for paragon in P x C",ylab="expression bias for charger in P x C")
abline(0,1)

## compares triad tpm between parents (r=0.95)
heatscatter(x=log2(PC_P_C$triad_tpm.x),y=log2(PC_P_C$triad_tpm.y),cor=T,main="",xlab="expression level (tpm) for paragon in P x C",ylab="expression level (tpm) for charger in P x C")
abline(0,1)

## is change in expression bias associated with in change in triad expression level? 
## weak correlation, so there seems to be compensation going on generally
heatscatter(x=log(PC_P_C$cv.x /PC_P_C$cv.y),y=log(PC_P_C$triad_tpm.x /PC_P_C$triad_tpm.y),cor=T,main="",xlab="log fold change in expression bias (cv)",ylab="log fold change in triad expression (tpm)")

dev.off()

pdf("PC parents biased triads.pdf",height=3,width=9)

## also apparent when looking at triads with high bias in one parent but not other. 
paragon_biased <- PC_P_C[PC_P_C$cv.x > 0.5 & PC_P_C$cv.y < 0.25, ]
p1 <- ggplot(data=paragon_biased,aes(x=log(cv.x /cv.y),y=log(triad_tpm.x /triad_tpm.y)))+
  geom_point() +
  xlab("log fold change in expression bias (cv)")+
  ylab("log fold change in triad expression (tpm)") +
  ggtitle("paragon_biased") +
  annotate("text", x = 1.5, y = -1.5, label = paste("n=",nrow(paragon_biased), ", R=",format(round(cor(log(paragon_biased$cv.x/paragon_biased$cv.y),log(paragon_biased$triad_tpm.x/paragon_biased$triad_tpm.y)), 2), nsmall = 3),sep=""))+
  geom_smooth(method=lm)
charger_biased <- PC_P_C[PC_P_C$cv.x < 0.25 & PC_P_C$cv.y > 0.5, ]
p2 <- ggplot(data=charger_biased,aes(x=log(cv.x /cv.y),y=log(triad_tpm.x /triad_tpm.y)))+
  geom_point() +
  xlab("log fold change in expression bias (cv)")+
  ylab("log fold change in triad expression (tpm)") +
  ggtitle("charger_biased") +
  annotate("text", x = -2.5, y = -0.5, label = paste("n=",nrow(charger_biased), ", R=",format(round(cor(log(charger_biased$cv.x/charger_biased$cv.y),log(charger_biased$triad_tpm.x/charger_biased$triad_tpm.y)), 2), nsmall = 3),sep=""))+
  geom_smooth(method=lm)

plot_grid(p1,p2,ncol=2)

dev.off()

## get all F5 lines only and add distance between parents
F6_PC <- HEB_PC[!HEB_PC$genotype %in% c("PC_P","PC_C"),]
F6_PC <- left_join(F6_PC,PC_P_C_sum,by="group_id")
F6_PC$ids <- paste(F6_PC$genotype,F6_PC$group_id,sep="_")

# only keep triads which are in both C and P
informative_group_ids <- HEB_PC %>%
  dplyr::select(group_id, genotype) %>%
  dplyr::filter(genotype == "PC_P" | genotype == "PC_C") %>% # only keep triads which are in W or P
  dplyr::group_by(group_id) %>%
  dplyr::summarise("in_both" = length(group_id)>1) %>%# assign "TRUE" to group IDs in both W and P
  dplyr::filter(in_both == TRUE)
F6_PC <- F6_PC %>%
  dplyr::filter(group_id %in% informative_group_ids$group_id)

## now calculate distance of each F5 line to P in PC
gids <- unique(F6_PC$genotype)
f6compareall <- ""
for (g in 1:length(gids)){
  f6compare <- inner_join(PC_P[c(2,8:10)],F6_PC[F6_PC$genotype %in% gids[g],c(2,8:10)],by="group_id")
  f6compare$distance_to_parent <- apply(f6compare[2:7], 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  f6compare <- f6compare[-c(2:7)]
  f6compare$genotype <- gids[g]
  f6compareall <- rbind(f6compareall,as.data.frame(f6compare))
}
f6compareall <- f6compareall[-c(1),]
f6compareall_distoP <- f6compareall
f6compareall_distoP$ids <- paste(f6compareall_distoP$genotype,f6compareall_distoP$group_id,sep="_")
f6compareall_distoP <- f6compareall_distoP[c(2,4)]
f6compareall_distoP$which_parent <- "To P"

## now calculate distance of each F5 line to C in PC
gids <- unique(F6_PC$genotype)
f6compareall <- ""
for (g in 1:length(gids)){
  f6compare <- inner_join(PC_C[c(2,8:10)],F6_PC[F6_PC$genotype %in% gids[g],c(2,8:10)],by="group_id")
  f6compare$distance_to_parent <- apply(f6compare[2:7], 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  f6compare <- f6compare[-c(2:7)]
  f6compare$genotype <- gids[g]
  f6compareall <- rbind(f6compareall,as.data.frame(f6compare))
}
f6compareall <- f6compareall[-c(1),]
f6compareall_distoC <- f6compareall
f6compareall_distoC$ids <- paste(f6compareall_distoC$genotype,f6compareall_distoC$group_id,sep="_")
f6compareall_distoC <- f6compareall_distoC[c(2,4)]
f6compareall_distoC$which_parent <- "To C"

## combine distances to compare
f6compareall_PC_distances <- inner_join(f6compareall_distoP,f6compareall_distoC,by="ids")
f6compareall_PC_distances$distance_to_parent.x <- as.numeric(f6compareall_PC_distances$distance_to_parent.x)
f6compareall_PC_distances$distance_to_parent.y <- as.numeric(f6compareall_PC_distances$distance_to_parent.y)

## add distance to parents to other F5 data
F6_PC1 <- left_join(F6_PC,f6compareall_distoP,by="ids")
F6_PC2 <- left_join(F6_PC,f6compareall_distoC,by="ids")
F6_PC3 <- rbind(F6_PC1,F6_PC2)
F6_PC3$distance_to_parent <- as.numeric(F6_PC3$distance_to_parent)
rm(F6_PC1)
rm(F6_PC2)

## no clear association between bias distance and expression bias, but maybe slight positive association (r = 0.22)
heatscatter(x=F6_PC3$distance_to_parent,y=F6_PC3$cv,main="",xlab = "bias distance",ylab="expression bias",cor=T)

## no clear association between bias distance and expression level, but maybe slight negative association (r = -0.21)
heatscatter(x=F6_PC3$distance_to_parent,y=log2(F6_PC3$triad_tpm),main="",xlab = "bias distance",ylab="expression level",cor=T)

## add parental cv and tpm as separate columns for comparison
F6_PC1 <- left_join(F6_PC,PC_P_cv_tpm,by="group_id")
F6_PC2 <- left_join(F6_PC,PC_C_cv_tpm,by="group_id")
F6_PCx <- rbind(F6_PC1,F6_PC2)
rm(F6_PC1)
rm(F6_PC2)

F6_PCwhole <- F6_PC3 %>% inner_join(F6_PCx)
write.table(F6_PCwhole,"F6_PCwhole",row.names = F,quote = F,sep = "\t")


jpeg("PC_comparisons.jpg",height=500,width=500)
par(mfrow = c(2, 2))

##  compare expression bias in parent and offspring - very strong correlation (0.77)
heatscatter(x=F6_PCx$cv_parent,y=F6_PCx$cv,main="",xlab = "expression bias in parents",ylab="expression bias in progeny", cor=T)
abline(0,1)

##  compare expression level in parent and offspring - very strong correlation (0.96)
heatscatter(x=log2(F6_PCx$tpm_parent),y=log2(F6_PCx$triad_tpm),main="",xlab = "expression level in parents",ylab="expression level in progeny", cor=T)
abline(0,1)

## comparing bias distance to parents, the two prongs represent full similarity to one parent and complete dissimilarity to other parent
heatscatter(x=F6_PC3$distance,y=F6_PC3$distance_to_parent,main="",xlab = "bias distance between parents",ylab="bias distance to parents")
#View(F6_PC3[F6_PC3$distance_to_parent < 0.1,])

## comparing bias distance to each parent, the two prongs represent full similarity to one parent and complete dissimilarity to other parent
heatscatter(x=f6compareall_PC_distances$distance_to_parent.x,y=f6compareall_PC_distances$distance_to_parent.y,main="",xlab = "distance to paragon in P x C",ylab="distance to charger in P x C")

dev.off()

## test for association between logFC in cv and logFC in tpm
cortest_pc <- F6_PCx %>% 
  dplyr::group_by(group_id) %>% 
  dplyr::summarise(pval = cor.test(x=log(cv/cv_parent),y=log(triad_tpm/tpm_parent))$p.value)
cortest_pc$adjustp <- p.adjust(cortest_pc$pval,method="holm")

## 1481 triads out of 2895 have significant association between expression bias and expression level
dim(cortest_pc)
dim(cortest_pc[cortest_pc$adjustp < 0.05,])
cortest_pc <- cortest_pc[order(cortest_pc$adjustp),]
tail(cortest_pc[cortest_pc$adjustp < 0.05,])

## example
ggplot(data=F6_PCx[F6_PCx$group_id %in% "5685",],aes(x=log(cv/cv_parent),y=log(triad_tpm/tpm_parent))) +
  geom_point() +
  xlab("log fold change in expression bias (cv) in F5 v parents") +
  ylab("log fold change in triad expression (tpm) in F5 v parents") +
  ggtitle("Example dosage sensitive triad") +
  geom_smooth(method=lm)

## now test for dosage compensation
## first using equivalence test to see if logfc in expression level in F5 and parents is significantly different from 0
eqitest_pc <- F6_PCx %>%
  dplyr::group_by(group_id) %>% 
  dplyr::summarise(pval = tost(log(triad_tpm/tpm_parent), 0, epsilon = 1, paired = F, var.equal = T, conf.level = 0.95)$tost.p.value)
eqitest_pc$adjustp <- p.adjust(eqitest_pc$pval,method="holm")
dim(eqitest_pc[eqitest_pc$adjustp < 0.05,])
#1081 dosage compensation

## now test if expression bias in F5 is significantly different from parents
onettest_pc <- F6_PCx %>%
  dplyr::group_by(group_id) %>% 
  dplyr::summarise(pval = t.test(log(cv/cv_parent),mu=0, alternative = "two.sided")$p.value)
onettest_pc$adjustp <- p.adjust(onettest_pc$pval,method="holm")
dim(onettest_pc[onettest_pc$adjustp < 0.05,])
#671 sig different from parents

## only keep triads with significance evidence of change in expression bias in F5 but equivalence of expression levels to parents
eqitest_pc2 <- eqitest_pc[eqitest_pc$group_id %in% onettest_pc[onettest_pc$adjustp < 0.05,]$group_id,]
dim(eqitest_pc2[eqitest_pc2$adjustp < 0.05,])
#301 kept

## only keep triads without evidence of association between change in expression level and change in expression bias
eqitest_pc2 <- eqitest_pc2[eqitest_pc2$adjustp < 0.05,][!eqitest_pc2[eqitest_pc2$adjustp < 0.05,]$group_id %in% cortest_pc[cortest_pc$adjustp < 0.05,]$group_id,]

#221 out of 2895 triads seem to have some evidence of dosage compensation
dim(eqitest_pc2)

ggplot(data=F6_PCx[F6_PCx$group_id %in% "18894",],aes(x=log(cv/cv_parent),y=log(triad_tpm/tpm_parent))) +
  geom_point() +
  xlab("log fold change in expression bias (cv) in F5 v parents") +
  ylab("log fold change in triad expression (tpm) in F5 v parents") +
  ggtitle("Example triad") +
  geom_smooth(method=lm)

# newdf_pc <- data.frame(Class = c("Dosage sensitive", "Dosage insensitive", "Other"), Proportion = c(1466/2895,203/2895,1226/2895))
# MAREK UPDATE - BASED ON THE NUMBERS FROM RUNNING ARUN'S SCRIPT

newdf_pc <- data.frame(Class = c("Dosage sensitive", "Dosage insensitive", "Other"), Proportion = c(1481/2951,221/2951,1249/2951))

ggplot(data=newdf_pc,aes(x=Class,y=Proportion)) +
  geom_col() +
  ggtitle("n = 2951 triads") +
  ylab("Proportion of triads") +
  xlab("Class of triads")

## comparing bias distance to each parent, the two prongs represent full similarity to one parent and complete dissimilarity to other parent
pdf("bias_distance_PC.pdf",height=4,width=4)
heatscatter(x=f6compareall_PC_distances$distance_to_parent.x,y=f6compareall_PC_distances$distance_to_parent.y,main="",xlab = "distance to paragon in P x C",ylab="distance to charger in P x C")
abline(0,1)
dev.off()

## are points in are different (in terms of expression bias distance) from both parents enrichment for some genotypes or triads?
f6compareall_PC_distances2 <- f6compareall_PC_distances %>%
  separate(ids,into=c("Cross","gt","group_id"))

## first study triads with extreme distance to both parents
extreme_distance_P_C <- table(f6compareall_PC_distances2[f6compareall_PC_distances2$distance_to_parent.x > 0.2 & f6compareall_PC_distances2$distance_to_parent.y > 0.2,]$group_id)
hist(extreme_distance_P_C,xlab="number of genotypes with extreme distance values",ylab="number of triads",main="")
length(extreme_distance_P_C)

## after isolating or removing triads with extreme distance to both parents, it still looks like many triads show increased distance from both parents, only for triads found in at least 15 genotypes
f6compareall_PC_distances_rm_with <- f6compareall_PC_distances2[f6compareall_PC_distances2$group_id %in% names(extreme_distance_P_C[extreme_distance_P_C > 15]),]
heatscatter(x=f6compareall_PC_distances_rm_with$distance_to_parent.x,y=f6compareall_PC_distances_rm_with$distance_to_parent.y,main="triads with extreme divergence from parents",xlab = "distance to paragon in P x C",ylab="distance to watkins in P x C")

## next study triads with extreme distance to paragon but not watkins
extreme_distance_P <- table(f6compareall_PC_distances2[f6compareall_PC_distances2$distance_to_parent.x > 0.2 & f6compareall_PC_distances2$distance_to_parent.y < 0.1,]$group_id)
hist(extreme_distance_P,xlab="number of genotypes with extreme distance values",ylab="number of triads",main="")
length(extreme_distance_P)

## after isolating or removing triads with extreme distance to paragon but not watkins, it still looks like many triads show increased distance from both parents
f6compareall_PC_distances_rm_with <- f6compareall_PC_distances2[f6compareall_PC_distances2$group_id %in% names(extreme_distance_P[extreme_distance_P > 15]),]
heatscatter(x=f6compareall_PC_distances_rm_with$distance_to_parent.x,y=f6compareall_PC_distances_rm_with$distance_to_parent.y,main="triads with extreme divergence from paragon only",xlab = "distance to paragon in P x C",ylab="distance to watkins in P x C")

## next study triads with extreme distance to watkins but not paragon
extreme_distance_C <- table(f6compareall_PC_distances2[f6compareall_PC_distances2$distance_to_parent.x < 0.1 & f6compareall_PC_distances2$distance_to_parent.y > 0.2,]$group_id)
hist(extreme_distance_C,xlab="number of genotypes with extreme distance values",ylab="number of triads",main="")
length(extreme_distance_C)

## after isolating or removing triads with extreme distance to paragon but not watkins, it still looks like many triads show increased distance from both parents
f6compareall_PC_distances_rm_with <- f6compareall_PC_distances2[f6compareall_PC_distances2$group_id %in% names(extreme_distance_C[extreme_distance_C > 15]),]
heatscatter(x=f6compareall_PC_distances_rm_with$distance_to_parent.x,y=f6compareall_PC_distances_rm_with$distance_to_parent.y,main="triads with extreme divergence from watkins only",xlab = "distance to paragon in P x C",ylab="distance to watkins in P x C")

## those that are taking extreme values in to just of one the parent are the same triads! 330 out of 513(P) or 499(C) triads

length(names(extreme_distance_P[extreme_distance_P > 15]))
length(names(extreme_distance_C[extreme_distance_C > 15]))
length(names(extreme_distance_P_C[extreme_distance_P_C > 15]))
sum(names(extreme_distance_P[extreme_distance_P > 15]) %in% names(extreme_distance_C[extreme_distance_C > 15]))

## but no sharing (only 5 triad shared!) between triads distance to only one parent and those distance to both parents 
sum(names(extreme_distance_P[extreme_distance_P > 15]) %in% names(extreme_distance_P_C[extreme_distance_P_C > 15]))
sum(names(extreme_distance_C[extreme_distance_C > 15]) %in% names(extreme_distance_P_C[extreme_distance_P_C > 15]))

## check distances of triads in parents
PC_P_C_sum_biased <- PC_P_C_sum
PC_P_C_sum_biased$divergence_in_F6 <- "other"
PC_P_C_sum_biased[PC_P_C_sum_biased$group_id %in% names(extreme_distance_P_C[extreme_distance_P_C > 15]),]$divergence_in_F6 <- "From both parents"
PC_P_C_sum_biased[PC_P_C_sum_biased$group_id %in% names(extreme_distance_P[extreme_distance_P > 15][which(names(extreme_distance_P[extreme_distance_P > 15]) %in% names(extreme_distance_C[extreme_distance_C > 15]))]),]$divergence_in_F6 <- "From one parent"
table(PC_P_C_sum_biased$divergence_in_F6)

ggplot(data=PC_P_C_sum_biased,aes(x=divergence_in_F6,y=distance))+
  geom_boxplot() +
  ylab("expression bias distance between parents")+
  xlab("triad bias divergence in F6")+ 
  geom_text(data=as.data.frame(table(PC_P_C_sum_biased$divergence_in_F6)),aes(x=Var1,y=max(PC_P_C_sum_biased$distance)+0.1,label=paste("n=", Freq)))

write.table(PC_P_C_sum_biased, file = "PC_P_C_sum_biased",quote = F,row.names = F,sep = "\t")


## now to PW specific analyses
## this part tests for triads with significant evidence of expression bias

## obtain gene expression data, note using all replicates, not just mean per genotype for stats
lmerdata <- TPM[grepl("PW_",TPM$genotype),]
lmerdata$sample <- gsub("PW_30.2.2","PW_30.3",lmerdata$sample) ## check whats going on with this sample with genewiz
lmerdata$sample <- gsub("\\.","-", lmerdata$sample)
lmerdata$sample <- gsub("-1-2$","-1",lmerdata$sample)
lmerdata$sample <- gsub("-2-2$","-2",lmerdata$sample)
lmerdata$sample <- gsub("-3-2$","-3",lmerdata$sample)
lmerdata$sample <- gsub("-4-2$","-4",lmerdata$sample)
lmerdata$replicate <- gsub(".*-","", lmerdata$sample )
lmerdata$group_id <- as.character(lmerdata$group_id)
rownames(lmerdata) <- paste(lmerdata$group_id,lmerdata$sample,sep=",")

## do a test for cv with genotype as fixed effect and replicate as random effect. 
## only for triads with at least 50 samples assayed
gentoypecvpvals_pw <- lmerdata[!lmerdata$group_id %in% c("16241","12652"),] %>%
  dplyr::group_by(group_id) %>%
  dplyr::filter(n()>50) %>%
  dplyr::summarise(pval = Anova(lmer(formula=cv~genotype+(1|replicate)))[,3])
## adjust pvalues using fdr to correct for multiple testing
gentoypecvpvals_pw$adjustP <- p.adjust(gentoypecvpvals_pw$pval,method = "BH")
gentoypecvpvals_pw$gt_sig <- T
gentoypecvpvals_pw[gentoypecvpvals_pw$adjustP > 0.05,]$gt_sig <- F
gentoypecvpvals_pw2 <- gentoypecvpvals_pw[-c(2:3)]
table(gentoypecvpvals_pw2$gt_sig)

## only analyze triads with significant evidence of expression bias variation between genotypes
HEB_PW <- HEB[grepl("PW_",HEB$genotype),]
HEB_PW <- HEB_PW[HEB_PW$group_id %in% gentoypecvpvals_pw[gentoypecvpvals_pw$adjustP < 0.0001,]$group_id,]

## parental triad expression estimates
PW_P <- HEB_PW %>%
  dplyr::filter(genotype == "PW_P")
PW_W <- HEB_PW %>%
  dplyr::filter(genotype == "PW_W")

PW_P_W <- inner_join(PW_P,PW_W,by="group_id")

PW_P_cv_tpm <- cbind(PW_P[c(2,11,12)],"To P")
colnames(PW_P_cv_tpm)[2:4] <- c("cv_parent","tpm_parent", "which_parent")
PW_W_cv_tpm <- cbind(PW_W[c(2,11,12)],"To W")
colnames(PW_W_cv_tpm)[2:4] <- c("cv_parent","tpm_parent", "which_parent")

## calculates distance matrix between parents
## this distance matrix is estimated by comparing first three columns to last three columns for each row separately
PW_all <- cbind(PW_P_W[8:10],PW_P_W[19:21])
PW_P_W$distance <- apply(PW_all, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
PW_P_W_sum <- cbind(PW_P_W[2],PW_P_W[ncol(PW_P_W)])

pdf("PW parents compare.pdf",height=3,width=9)
par(mfrow = c(1, 3))

## compares HEB between parents (r=0.75)
heatscatter(x=PW_P_W$cv.x,y=PW_P_W$cv.y,cor=T,main="",xlab="expression bias for paragon in P x W",ylab="expression bias for watkins in P x W")
abline(0,1)

## compares triad tpm between parents (r=0.95)
heatscatter(x=log2(PW_P_W$triad_tpm.x),y=log2(PW_P_W$triad_tpm.y),cor=T,main="",xlab="expression level (tpm) for paragon in P x W",ylab="expression level (tpm) for watkins in P x W")
abline(0,1)

## is change in expression bias associated with in change in triad expression level? 
## weak correlation, so there seems to be compensation going on generally
heatscatter(x=log(PW_P_W$cv.x /PW_P_W$cv.y),y=log(PW_P_W$triad_tpm.x /PW_P_W$triad_tpm.y),cor=T,main="",xlab="log fold change in expression bias (cv)",ylab="log fold change in triad expression (tpm)")

dev.off()

pdf("PW parents biased triads.pdf",height=3,width=9)

## also apparent when looking at triads with high bias in one parent but not other. 
paragon_biased <- PW_P_W[PW_P_W$cv.x > 0.5 & PW_P_W$cv.y < 0.25, ]
p1 <- ggplot(data=paragon_biased,aes(x=log(cv.x /cv.y),y=log(triad_tpm.x /triad_tpm.y)))+
  geom_point() +
  xlab("log fold change in expression bias (cv)")+
  ylab("log fold change in triad expression (tpm)") +
  ggtitle("paragon_biased") +
  annotate("text", x = 1.5, y = -1.5, label = paste("n=",nrow(paragon_biased), ", R=",format(round(cor(log(paragon_biased$cv.x/paragon_biased$cv.y),log(paragon_biased$triad_tpm.x/paragon_biased$triad_tpm.y)), 2), nsmall = 3),sep=""))+
  geom_smooth(method=lm)
watkins_biased <- PW_P_W[PW_P_W$cv.x < 0.25 & PW_P_W$cv.y > 0.5, ]
p2 <- ggplot(data=watkins_biased,aes(x=log(cv.x /cv.y),y=log(triad_tpm.x /triad_tpm.y)))+
  geom_point() +
  xlab("log fold change in expression bias (cv)")+
  ylab("log fold change in triad expression (tpm)") +
  ggtitle("watkins_biased") +
  annotate("text", x = -2.5, y = -0.5, label = paste("n=",nrow(watkins_biased), ", R=",format(round(cor(log(watkins_biased$cv.x/watkins_biased$cv.y),log(watkins_biased$triad_tpm.x/watkins_biased$triad_tpm.y)), 2), nsmall = 3),sep=""))+
  geom_smooth(method=lm)

plot_grid(p1,p2,ncol=2)

dev.off()

## get all F5 lines only and add distance between parents
F6_PW <- HEB_PW[!HEB_PW$genotype %in% c("PW_P","PW_W"),]
F6_PW <- left_join(F6_PW,PW_P_W_sum,by="group_id")
F6_PW$ids <- paste(F6_PW$genotype,F6_PW$group_id,sep="_")

# only keep triads which are in both W and P
informative_group_ids <- HEB_PW %>%
  dplyr::select(group_id, genotype) %>%
  dplyr::filter(genotype == "PW_P" | genotype == "PW_W") %>% # only keep triads which are in W or P
  dplyr::group_by(group_id) %>%
  dplyr::summarise("in_both" = length(group_id)>1) %>%# assign "TRUE" to group IDs in both W and P
  dplyr::filter(in_both == TRUE)
F6_PW <- F6_PW %>%
  dplyr::filter(group_id %in% informative_group_ids$group_id)

## now calculate distance of each F5 line to P in PW
gids <- unique(F6_PW$genotype)
f6compareall <- ""
for (g in 1:length(gids)){
  f6compare <- inner_join(PW_P[c(2,8:10)],F6_PW[F6_PW$genotype %in% gids[g],c(2,8:10)],by="group_id")
  f6compare$distance_to_parent <- apply(f6compare[2:7], 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  f6compare <- f6compare[-c(2:7)]
  f6compare$genotype <- gids[g]
  f6compareall <- rbind(f6compareall,as.data.frame(f6compare))
}
f6compareall <- f6compareall[-c(1),]
f6compareall_distoP <- f6compareall
f6compareall_distoP$ids <- paste(f6compareall_distoP$genotype,f6compareall_distoP$group_id,sep="_")
f6compareall_distoP <- f6compareall_distoP[c(2,4)]
f6compareall_distoP$which_parent <- "To P"

## now calculate distance of each F5 line to W in PW
gids <- unique(F6_PW$genotype)
f6compareall <- ""
for (g in 1:length(gids)){
  f6compare <- inner_join(PW_W[c(2,8:10)],F6_PW[F6_PW$genotype %in% gids[g],c(2,8:10)],by="group_id")
  f6compare$distance_to_parent <- apply(f6compare[2:7], 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
  f6compare <- f6compare[-c(2:7)]
  f6compare$genotype <- gids[g]
  f6compareall <- rbind(f6compareall,as.data.frame(f6compare))
}
f6compareall <- f6compareall[-c(1),]
f6compareall_distoW <- f6compareall
f6compareall_distoW$ids <- paste(f6compareall_distoW$genotype,f6compareall_distoW$group_id,sep="_")
f6compareall_distoW <- f6compareall_distoW[c(2,4)]
f6compareall_distoW$which_parent <- "To W"

## combine distances to compare
f6compareall_PW_distances <- inner_join(f6compareall_distoP,f6compareall_distoW,by="ids")
f6compareall_PW_distances$distance_to_parent.x <- as.numeric(f6compareall_PW_distances$distance_to_parent.x)
f6compareall_PW_distances$distance_to_parent.y <- as.numeric(f6compareall_PW_distances$distance_to_parent.y)

## add distance to parents to other F5 data
F6_PW1 <- left_join(F6_PW,f6compareall_distoP,by="ids")
F6_PW2 <- left_join(F6_PW,f6compareall_distoW,by="ids")
F6_PW3 <- rbind(F6_PW1,F6_PW2)
F6_PW3$distance_to_parent <- as.numeric(F6_PW3$distance_to_parent)
rm(F6_PW1)
rm(F6_PW2)

## no association between bias distance and expression bias, but maybe slight positive association (r = 0.22)
#heatscatter(x=F6_PW3$distance_to_parent,y=F6_PW3$cv,main="",xlab = "bias distance",ylab="expression bias",cor=T)

## no association between bias distance and expression level, but maybe slight negative association (r = -0.21)
#heatscatter(x=F6_PW3$distance_to_parent,y=log2(F6_PW3$triad_tpm),main="",xlab = "bias distance",ylab="expression level",cor=T)

## add parental cv and tpm as separate columns for comparison
F6_PW1 <- left_join(F6_PW,PW_P_cv_tpm,by="group_id")
F6_PW2 <- left_join(F6_PW,PW_W_cv_tpm,by="group_id")
F6_PWx <- rbind(F6_PW1,F6_PW2)
rm(F6_PW1)
rm(F6_PW2)

F6_PWwhole <- F6_PW3 %>% inner_join(F6_PWx)
write.table(F6_PWwhole,"F6_PWwhole",row.names = F,quote = F,sep = "\t")


jpeg("PW_comparisons.jpg",height=500,width=500)
par(mfrow = c(2, 2))

##  compare expression bias in parent and offspring - very strong correlation (0.77)
heatscatter(x=F6_PWx$cv_parent,y=F6_PWx$cv,main="",xlab = "expression bias in parents",ylab="expression bias in progeny", cor=T)
abline(0,1)

##  compare expression level in parent and offspring - very strong correlation (0.96)
heatscatter(x=log2(F6_PWx$tpm_parent),y=log2(F6_PWx$triad_tpm),main="",xlab = "expression level in parents",ylab="expression level in progeny", cor=T)
abline(0,1)

## comparing bias distance to parents, the two prongs represent full similarity to one parent and complete dissimilarity to other parent
heatscatter(x=F6_PW3$distance,y=F6_PW3$distance_to_parent,main="",xlab = "bias distance between parents",ylab="bias distance to parents")
#View(F6_PW3[F6_PW3$distance_to_parent < 0.1,])

## comparing bias distance to each parent, the two prongs represent full similarity to one parent and complete dissimilarity to other parent
heatscatter(x=f6compareall_PW_distances$distance_to_parent.x,y=f6compareall_PW_distances$distance_to_parent.y,main="",xlab = "distance to paragon in P x W",ylab="distance to watkins in P x W")

dev.off()

pdf("bias_distance_PW.pdf",height = 4,width = 4)
heatscatter(x=f6compareall_PW_distances$distance_to_parent.x,y=f6compareall_PW_distances$distance_to_parent.y,main="",xlab = "distance to paragon in P x W",ylab="distance to watkins in P x W")
dev.off()

## test for association between logFC in cv and logFC in tpm
cortest_pw <- F6_PWx %>% 
  dplyr::group_by(group_id) %>% 
  dplyr::summarise(pval = cor.test(x=log(cv/cv_parent),y=log(triad_tpm/tpm_parent))$p.value)
cortest_pw$adjustp <- p.adjust(cortest_pw$pval,method="holm")

## 1112 triads out of 2663 have significant association between expression bias and expression level

dim(cortest_pw)
dim(cortest_pw[cortest_pw$adjustp < 0.05,])
cortest_pw <- cortest_pw[order(cortest_pw$adjustp),]
tail(cortest_pw[cortest_pw$adjustp < 0.05,])

## example
ggplot(data=F6_PWx[F6_PWx$group_id %in% "5261",],aes(x=log(cv/cv_parent),y=log(triad_tpm/tpm_parent))) +
  geom_point() +
  xlab("log fold change in expression\nbias (cv) in F5 v parents") +
  ylab("log fold change in triad\nexpression (tpm) in F5 v parents") +
  ggtitle("Example triad") +
  geom_smooth(method=lm)

pdf("dosage_sensitive_triad.pdf",height=3.4,width=3.4)
ggplot(data=F6_PCx[F6_PCx$group_id %in% "5261",],aes(x=log(cv/cv_parent),y=log(triad_tpm/tpm_parent))) +
  geom_point() +
  xlab("log fold change in expression\nbias (cv) in F5 v parents") +
  ylab("log fold change in triad\nexpression (tpm) in F5 v parents") +
  ggtitle("Example dosage sensitive triad") +
  geom_smooth(method=lm)+
  theme_bw()
dev.off()

## now test for dosage compensation
## first using equivalence test to see if logfc in expression level in F5 and parents is significantly different from 0
eqitest_pw <- F6_PWx %>%
  dplyr::group_by(group_id) %>% 
  dplyr::summarise(pval = tost(log(triad_tpm/tpm_parent), 0, epsilon = 1, paired = F, var.equal = T, conf.level = 0.95)$tost.p.value)
eqitest_pw$adjustp <- p.adjust(eqitest_pw$pval,method="holm")
dim(eqitest_pw[eqitest_pw$adjustp < 0.05,])
#487

## now test if expression bias in F5 is significantly different from parents
onettest_pw <- F6_PWx %>%
  dplyr::group_by(group_id) %>% 
  dplyr::summarise(pval = t.test(log(cv/cv_parent),mu=0, alternative = "two.sided")$p.value)
onettest_pw$adjustp <- p.adjust(onettest_pw$pval,method="holm")
dim(onettest_pw[onettest_pw$adjustp < 0.05,])
#1284

## only keep triads with significance evidence of change in expression bias in F5 but equivalence of expression levels to parents
eqitest_pw2 <- eqitest_pw[eqitest_pw$group_id %in% onettest_pw[onettest_pw$adjustp < 0.05,]$group_id,]
dim(eqitest_pw2[eqitest_pw2$adjustp < 0.05,])
#263

## only keep triads without evidence of association between change in expression level and change in expression bias
eqitest_pw2 <- eqitest_pw2[eqitest_pw2$adjustp < 0.05,][!eqitest_pw2[eqitest_pw2$adjustp < 0.05,]$group_id %in% cortest_pw[cortest_pw$adjustp < 0.05,]$group_id,]

#172 out of 2663 triads seem to have some evidence of dosage compensation

dim(eqitest_pw2)

ggplot(data=F6_PWx[F6_PWx$group_id %in% "6265",],aes(x=log(cv/cv_parent,),y=log(triad_tpm/tpm_parent))) +
  geom_point() +
  xlab("log fold change in expression\nbias (cv) in F5 v parents") +
  ylab("log fold change in triad\nexpression (tpm) in F5 v parents") +
  ggtitle("Example dosage insensitive triad") +
  geom_smooth(method=lm)

pdf("dosage_insensitive_triad.pdf",height=3.4,width=3.4)
ggplot(data=F6_PWx[F6_PWx$group_id %in% "6265",],aes(x=log(cv/cv_parent,),y=log(triad_tpm/tpm_parent))) +
  geom_point() +
  xlab("log fold change in expression\nbias (cv) in F5 v parents") +
  ylab("log fold change in triad\nexpression (tpm) in F5 v parents") +
  ggtitle("Example dosage insensitive triad") +
  geom_smooth(method=lm)+
  theme_bw()
dev.off()

newdf_pw <- data.frame(Class = c("Dosage\nsensitive", "Dosage\ninsensitive", "Other"), Proportion = c(1112/2663,172/2663,1379/2663))

pdf("proportion_triads.pdf",height=2.5,width=2.5)
ggplot(data=newdf_pw,aes(x=Class,y=Proportion)) +
  geom_col() +
  ggtitle("n = 2663 triads") +
  ylab("Proportion of triads") +
  xlab("")+
  theme_bw()
dev.off()

## comparing bias distance to each parent, the two prongs represent full similarity to one parent and complete dissimilarity to other parent
pdf("bias_distance_PW.pdf",height=4,width=4)
heatscatter(x=f6compareall_PW_distances$distance_to_parent.x,y=f6compareall_PW_distances$distance_to_parent.y,main="",xlab = "distance to paragon in P x W",ylab="distance to watkins in P x W")
abline(0,1)
dev.off()

## are points in are different (in terms of expression bias distance) from both parents enrichment for some genotypes or triads?
f6compareall_PW_distances2 <- f6compareall_PW_distances %>%
  separate(ids,into=c("Cross","gt","group_id"))

## first study triads with extreme distance to both parents
extreme_distance_P_W <- table(f6compareall_PW_distances2[f6compareall_PW_distances2$distance_to_parent.x > 0.2 & f6compareall_PW_distances2$distance_to_parent.y > 0.2,]$group_id)
hist(extreme_distance_P_W,xlab="number of genotypes with extreme distance values",ylab="number of triads",main="")
length(extreme_distance_P_W)
f6compareall_PW_distances_rm_with <- f6compareall_PW_distances2[f6compareall_PW_distances2$group_id %in% names(extreme_distance_P_W[extreme_distance_P_W > 15]),]
heatscatter(x=f6compareall_PW_distances_rm_with$distance_to_parent.x,y=f6compareall_PW_distances_rm_with$distance_to_parent.y,main="triads with extreme divergence from both parents",xlab = "distance to paragon in P x W",ylab="distance to watkins in P x W")

pdf("highly_diverged_bias_distance_PW.pdf",height=4,width=3.5)
heatscatter(x=f6compareall_PW_distances_rm_with$distance_to_parent.x,y=f6compareall_PW_distances_rm_with$distance_to_parent.y,main="",xlab = "Bias distance to parent 1",ylab="Bias distance to parent 2")
abline(0,1)
dev.off()

## What do triads distant from both parents look for example triads
names(extreme_distance_P_W[extreme_distance_P_W > 15])
heatscatter(x=f6compareall_PW_distances_rm_with[f6compareall_PW_distances_rm_with$group_id %in% "10170",]$distance_to_parent.x,y=f6compareall_PW_distances_rm_with[f6compareall_PW_distances_rm_with$group_id %in% "10170",]$distance_to_parent.y,main="triads with extreme divergence from watkins only",xlab = "distance to paragon in P x W",ylab="distance to watkins in P x W")
HEB_PW_plot <- as.data.frame(HEB_PW[HEB_PW$group_id %in% 10170,])
rownames(HEB_PW_plot) <- HEB_PW_plot$genotype
pheatmap(t(HEB_PW_plot[,8:10]))

pdf(file="divergent_both_parents.pdf",width=5.2,height=3)
pheatmap(t(HEB_PW_plot[,8:10]),cluster_rows = F,labels_col = c(rep("",49),rep("P",2)))
dev.off()

heatscatter(x=f6compareall_PW_distances_rm_with[f6compareall_PW_distances_rm_with$group_id %in% "10256",]$distance_to_parent.x,y=f6compareall_PW_distances_rm_with[f6compareall_PW_distances_rm_with$group_id %in% "10256",]$distance_to_parent.y,main="triads with extreme divergence from watkins only",xlab = "distance to paragon in P x W",ylab="distance to watkins in P x W")
HEB_PW_plot <- as.data.frame(HEB_PW[HEB_PW$group_id %in% 10256,])
rownames(HEB_PW_plot) <- HEB_PW_plot$genotype
pheatmap(t(HEB_PW_plot[,8:10]))

## next study triads with extreme distance to paragon but not watkins
extreme_distance_P <- table(f6compareall_PW_distances2[f6compareall_PW_distances2$distance_to_parent.x > 0.2 & f6compareall_PW_distances2$distance_to_parent.y < 0.1,]$group_id)
hist(extreme_distance_P,xlab="number of genotypes with extreme distance values",ylab="number of triads",main="")
length(extreme_distance_P)

## after isolating or removing triads with extreme distance to paragon but not watkins, it still looks like many triads show increased distance from both parents
f6compareall_PW_distances_rm_with <- f6compareall_PW_distances2[f6compareall_PW_distances2$group_id %in% names(extreme_distance_P[extreme_distance_P > 15]),]
heatscatter(x=f6compareall_PW_distances_rm_with$distance_to_parent.x,y=f6compareall_PW_distances_rm_with$distance_to_parent.y,main="triads with extreme divergence from paragon only",xlab = "distance to paragon in P x W",ylab="distance to watkins in P x W")

## next study triads with extreme distance to watkins but not paragon
extreme_distance_W <- table(f6compareall_PW_distances2[f6compareall_PW_distances2$distance_to_parent.x < 0.1 & f6compareall_PW_distances2$distance_to_parent.y > 0.2,]$group_id)
hist(extreme_distance_W,xlab="number of genotypes with extreme distance values",ylab="number of triads",main="")
length(extreme_distance_W)

## after isolating or removing triads with extreme distance to watkins but not paragon, it still looks like many triads show increased distance from both parents
f6compareall_PW_distances_rm_with <- f6compareall_PW_distances2[f6compareall_PW_distances2$group_id %in% names(extreme_distance_W[extreme_distance_W > 15]),]
heatscatter(x=f6compareall_PW_distances_rm_with$distance_to_parent.x,y=f6compareall_PW_distances_rm_with$distance_to_parent.y,main="triads with extreme divergence from watkins only",xlab = "distance to paragon in P x W",ylab="distance to watkins in P x W")

pdf("highly_diverged_bias_watkins_PW.pdf",height=4,width=3.5)
heatscatter(x=f6compareall_PW_distances_rm_with$distance_to_parent.x,y=f6compareall_PW_distances_rm_with$distance_to_parent.y,main="",xlab = "Bias distance to parent 1",ylab="Bias distance to parent 2")
dev.off()

## those that are taking extreme values in to just of one the parent are the same triads! 153 out of 227 or 219 triads
## those that are taking extreme values in to just of one the parent are the same triads! 155 out of 220 (P) or 229 (W) triads
length(names(extreme_distance_P[extreme_distance_P > 15]))
length(names(extreme_distance_W[extreme_distance_W > 15]))
length(names(extreme_distance_P_W[extreme_distance_P_W > 15]))
sum(names(extreme_distance_P[extreme_distance_P > 15]) %in% names(extreme_distance_W[extreme_distance_W > 15]))

## What does parent matching (distant from one parent but close to the other) look for example triads
names(extreme_distance_P[extreme_distance_P > 15]) [names(extreme_distance_P[extreme_distance_P > 15]) %in% names(extreme_distance_W[extreme_distance_W > 15])]
heatscatter(x=f6compareall_PW_distances_rm_with[f6compareall_PW_distances_rm_with$group_id %in% "10117",]$distance_to_parent.x,y=f6compareall_PW_distances_rm_with[f6compareall_PW_distances_rm_with$group_id %in% "10117",]$distance_to_parent.y,main="triads with extreme divergence from watkins only",xlab = "distance to paragon in P x W",ylab="distance to watkins in P x W")
HEB_PW_plot <- as.data.frame(HEB_PW[HEB_PW$group_id %in% 10117,])
rownames(HEB_PW_plot) <- HEB_PW_plot$genotype
pheatmap(t(HEB_PW_plot[,8:10]))
heatscatter(x=f6compareall_PW_distances_rm_with[f6compareall_PW_distances_rm_with$group_id %in% "10175",]$distance_to_parent.x,y=f6compareall_PW_distances_rm_with[f6compareall_PW_distances_rm_with$group_id %in% "10175",]$distance_to_parent.y,main="triads with extreme divergence from watkins only",xlab = "distance to paragon in P x W",ylab="distance to watkins in P x W")
HEB_PW_plot <- as.data.frame(HEB_PW[HEB_PW$group_id %in% 10175,])
rownames(HEB_PW_plot) <- HEB_PW_plot$genotype
pheatmap(t(HEB_PW_plot[,8:10]))

pdf(file="divergent_one_parent.pdf",width=5.2,height=3)
pheatmap(t(HEB_PW_plot[,8:10]),cluster_rows = F,labels_col = c(rep("",50),rep("P",2)))
dev.off()

## but no sharing (only 1 triad shared!) between triads distance to only one parent and those distance to both parents 
sum(names(extreme_distance_P[extreme_distance_P > 15]) %in% names(extreme_distance_P_W[extreme_distance_P_W > 15]))
sum(names(extreme_distance_W[extreme_distance_W > 15]) %in% names(extreme_distance_P_W[extreme_distance_P_W > 15]))

## check distances of triads in parents
PW_P_W_sum_biased <- PW_P_W_sum
PW_P_W_sum_biased$divergence_in_F6 <- "other"
PW_P_W_sum_biased[PW_P_W_sum_biased$group_id %in% names(extreme_distance_P_W[extreme_distance_P_W > 15]),]$divergence_in_F6 <- "From both parents"
PW_P_W_sum_biased[PW_P_W_sum_biased$group_id %in% names(extreme_distance_P[extreme_distance_P > 15][which(names(extreme_distance_P[extreme_distance_P > 15]) %in% names(extreme_distance_W[extreme_distance_W > 15]))]),]$divergence_in_F6 <- "From one parent"
table(PW_P_W_sum_biased$divergence_in_F6)
# PW_P_W_sum_biased$divergence_in_F6 <- gsub("From both parents","From both parents",PW_P_W_sum_biased$divergence_in_F6)
# PW_P_W_sum_biased$divergence_in_F6 <- gsub("From one parent","From one parent",PW_P_W_sum_biased$divergence_in_F6)

pdf(file="parent_distances.pdf",width=2.8,height=2.8)
ggplot(data=PW_P_W_sum_biased,aes(x=divergence_in_F6,y=distance))+
  geom_boxplot() +
  ylab("Expression bias distance\nbetween parents")+
  xlab("Triad bias divergence in F6")+ 
  theme_bw() +
  geom_text(data=as.data.frame(table(PW_P_W_sum_biased$divergence_in_F6)),aes(x=Var1,y=max(PW_P_W_sum_biased$distance)+0.1,label=paste("n=", Freq)))
dev.off()

write.table(PW_P_W_sum_biased,file = "PW_P_W_sum_biased",row.names = F,quote = F,sep = "\t")
