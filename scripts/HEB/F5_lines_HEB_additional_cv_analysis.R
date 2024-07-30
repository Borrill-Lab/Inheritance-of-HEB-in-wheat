# Marek Glombik
# This script revises the analysis based on CV to show how it improves HEB analysis instead of using categories

library(tidyverse)
library(dplyr)
library(pracma)
library(ggplot2)
library(ggtern)
library(svglite)
library(ggpattern)
library(ggpubr)

#Check how HEB varies amongst lines
bias_categ_f6 <- read.csv('bias_category_all_samples_inc_orig_expr.csv',header = T)

bias_categ_f6$f6line <- substr(bias_categ_f6$sample,1,5)

#Quick info about how many triads per sample were detected + average in each parent

head(bias_categ_f6[bias_categ_f6$sample=='PC.1.1',])
quicksum <- bias_categ_f6 %>%
  group_by(sample) %>%
  summarise(triads_expressed=n())
min(quicksum$triads_expressed) #11,215
max(quicksum$triads_expressed) #14,331
subset_quicksum <- quicksum[grepl("PW.W|PW.P|PC.P|PC.C", quicksum$sample), ]
mean(subset_quicksum$triads_expressed)

#Filter out the parents, because we do not want them in this summary

bias_categ_f6_nopar <- bias_categ_f6[bias_categ_f6$f6line!='PC.P.' & bias_categ_f6$f6line!='PC.C.' & 
                                       bias_categ_f6$f6line!='PW.P.' & bias_categ_f6$f6line!='PW.W.',]


bias_categ_f6_sum <- bias_categ_f6_nopar %>%
  group_by(group_id,f6line) %>%
  summarise(distinct_in_triad=n_distinct(name_mins))

checksum_distinct <- bias_categ_f6_sum[bias_categ_f6_sum$group_id==2,]
checksum_bias <- bias_categ_f6[bias_categ_f6$group_id==2,]


bias_categ_f6_sum$cross <- substr(bias_categ_f6_sum$f6line,1,2)

total_samples_triad <- bias_categ_f6_sum %>%
  group_by(cross,group_id) %>%
  summarise(total_samples=n())

bias_categ_f6_sum_sum <- bias_categ_f6_sum %>%
  group_by(cross,group_id,distinct_in_triad) %>%
  summarise(total_distinct=n()) %>%
  inner_join(total_samples_triad)

bias_categ_f6_sum_sum$perc <- bias_categ_f6_sum_sum$total_distinct/bias_categ_f6_sum_sum$total_samples*100

#How many triads we have in total info for
dim(unique(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$cross=='PC',2]))
# 16411
dim(unique(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$cross=='PW',2]))
# 16241

#How many are represented by at least 50 samples? 

dim(unique(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$cross=='PC' & bias_categ_f6_sum_sum$total_samples==50,2]))
# 12717
dim(unique(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$cross=='PW' & bias_categ_f6_sum_sum$total_samples==50,2]))
# 12598

# And present in both parents?
dim(unique(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$cross=='PC' & bias_categ_f6_sum_sum$group_id %in% triadsPCpar$group_id,c('group_id')]))



dim(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$distinct_in_triad>1 & bias_categ_f6_sum_sum$cross=='PC' & bias_categ_f6_sum_sum$total_samples==50,])

dim(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$distinct_in_triad==1 & bias_categ_f6_sum_sum$cross=='PC' & bias_categ_f6_sum_sum$total_samples==50,])

#How many are labelled in 1 category in 100% of cases?
dim(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$distinct_in_triad==1 & bias_categ_f6_sum_sum$cross=='PC' & 
                            bias_categ_f6_sum_sum$total_samples==50 & bias_categ_f6_sum_sum$total_distinct==50,])
# 2929

dim(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$distinct_in_triad==1 & bias_categ_f6_sum_sum$cross=='PW' & 
                            bias_categ_f6_sum_sum$total_samples==50 & bias_categ_f6_sum_sum$total_distinct==50,])
# 3179

#How many are >=90 % = 45 out of 50 samples?
dim(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$distinct_in_triad==1 & bias_categ_f6_sum_sum$cross=='PC' & 
                            bias_categ_f6_sum_sum$total_samples==50 & bias_categ_f6_sum_sum$perc>=90,])
# 5519

dim(bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$distinct_in_triad==1 & bias_categ_f6_sum_sum$cross=='PW' & 
                            bias_categ_f6_sum_sum$total_samples==50 & bias_categ_f6_sum_sum$perc>=90,])
# 5731


# Try the line graph as Phlippa suggested for it

over_50_samp_f6_distinct <- bias_categ_f6_sum_sum[bias_categ_f6_sum_sum$total_samples==50,]
over_50_samp_f6_distinct <- over_50_samp_f6_distinct %>%
  arrange(cross,perc,distinct_in_triad) %>%
  group_by(perc,cross,distinct_in_triad) %>%
  summarise(countperc=n())

over_50_samp_f6_distinct$distinct_in_triad <- as.character(over_50_samp_f6_distinct$distinct_in_triad)
colnames(over_50_samp_f6_distinct)[3] <- 'Categories' 
ggplot(over_50_samp_f6_distinct,aes(perc,countperc,color=cross,linetype=Categories)) +
  geom_line(linewidth=1,alpha=0.8,position = position_dodge(width = 8)) +
  scale_linetype_manual(values = c("solid","dashed", "dotted")) +
  scale_color_manual(values = c("black", "#009128ff"),labels=c('PxC','PxW'),name='Cross') +
  theme_classic() +
  xlab(expression('Percentage of F'[5]~'lines with biological replicates in 1, 2, or 3 categories')) +
  ylab('Number of triads') +
  annotate('text',x=70,y=3200,label=paste0('PxW Total triads = 12,598\n100 % in 1 category = 3,179 (25 %)'),color='#009128ff',size=7) +
  annotate('text',x=65,y=2800,label=paste0('PxC Total triads = 12,717\n100 % in 1 category = 2,929 (23 %)'),color='black',size=7) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size=22),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 22))
  

ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_line_cumulative.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_line_cumulative.png',device = 'png',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

### Supp figure - plot it as a triangular graph showing the density of triads based on the multiple categories

trg_bias_categ <- bias_categ_f6_nopar %>%
  inner_join(bias_categ_f6_sum)

trg_bias_categPC <- trg_bias_categ[trg_bias_categ$cross=='PC',]
trg_bias_categPC_dist1PC10 <- trg_bias_categPC[trg_bias_categPC$distinct_in_triad==1 & trg_bias_categPC$f6line=='PC.35',]
trg_bias_categPC_dist3PC10 <- trg_bias_categPC[trg_bias_categPC$distinct_in_triad==3 & trg_bias_categPC$f6line=='PC.35',]
trg_bias_categPC_dist2PC10 <- trg_bias_categPC[trg_bias_categPC$distinct_in_triad==2 & trg_bias_categPC$f6line=='PC.35',]

ggtern(trg_bias_categPC_dist2PC10,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='black') +
  scale_fill_gradient2(low = 'white',high = 'black') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_2catPC.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_2catPC.png',device = 'png',height = 8,width = 12,dpi = 800,units = 'in',plot = last_plot())
ggtern(trg_bias_categPC_dist3PC10,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='black') +
  scale_fill_gradient2(low = 'white',high = 'black') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_3catPC.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_3catPC.png',device = 'png',height = 8,width = 12,dpi = 800,units = 'in',plot = last_plot())
ggtern(trg_bias_categPC_dist1PC10,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='black') +
  scale_fill_gradient2(low = 'white',high = 'black') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_1catPC.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_1catPC.png',device = 'png',height = 8,width = 12,dpi = 800,units = 'in',plot = last_plot())


trg_bias_categPW <- trg_bias_categ[trg_bias_categ$cross=='PW',]
trg_bias_categPW_dist1PW10 <- trg_bias_categPW[trg_bias_categPW$distinct_in_triad==1 & trg_bias_categPW$f6line=='PW.17',]
trg_bias_categPW_dist3PW10 <- trg_bias_categPW[trg_bias_categPW$distinct_in_triad==3 & trg_bias_categPW$f6line=='PW.17',]
trg_bias_categPW_dist2PW10 <- trg_bias_categPW[trg_bias_categPW$distinct_in_triad==2 & trg_bias_categPW$f6line=='PW.17',]

ggtern(trg_bias_categPW_dist2PW10,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='#009128ff') +
  scale_fill_gradient2(low = 'white',high = '#009128ff') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_2catPW.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_2catPW.png',device = 'png',height = 8,width = 12,dpi = 800,units = 'in',plot = last_plot())
ggtern(trg_bias_categPW_dist3PW10,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='#009128ff') +
  scale_fill_gradient2(low = 'white',high = '#009128ff') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_3catPW.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_3catPW.png',device = 'png',height = 8,width = 12,dpi = 800,units = 'in',plot = last_plot())
ggtern(trg_bias_categPW_dist1PW10,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='#009128ff') +
  scale_fill_gradient2(low = 'white',high = '#009128ff') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_1catPW.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_in_reps_ternary_1catPW.png',device = 'png',height = 8,width = 12,dpi = 800,units = 'in',plot = last_plot())


#####
# Next is to select a sample triad with reps in different categories to represent them on the schematic figure
# Best some with a low area spread but crossing the line

bias2distinct <- bias_categ_f6_sum[bias_categ_f6_sum$distinct_in_triad==2,]

schemeCV_PC <- bias_categ_f6[bias_categ_f6$f6line=='PW.11' & bias_categ_f6$group_id=='8571',]
schemeCV_PC$cv <- apply(schemeCV_PC[4:6],1,sd)/apply(schemeCV_PC[4:6],1,mean)
mean(schemeCV_PC$cv)


#Now plot it in a triangle

ggtern(schemeCV_PC,aes(A_tpm,D_tpm,B_tpm)) +
  geom_Tline(Tintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_Rline(Rintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_Lline(Lintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_mask() +
  geom_point(colour="black",size=3,alpha=0.8,fill='white',pch=21) +
  labs(x="A",y="D",z="B") +
  theme_gray() +
  theme(axis.title = element_text(size = 15),axis.text = element_text(size = 10))
ggsave('/Users/glombik/work/obj1_reanalysis/revised/HEB_variation_triangle.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())


#################
# Reuse Arun's code from F6_lines_heb.R to make additional plots - do the cv~f6line and plot summary for both populations


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
## only analyze triads with significant evidence of expression bias variation between genotypes
HEB_PC <- HEB[grepl("PC_",HEB$genotype),]
HEB_PC <- HEB_PC[HEB_PC$group_id %in% gentoypecvpvals_pc[gentoypecvpvals_pc$adjustP < 0.0001,]$group_id,]

gentoypecvpvals_pc_0001 <- gentoypecvpvals_pc[gentoypecvpvals_pc$adjustP<0.0001,]

# FOr the plot I need those that are in both parents (thats how I will get 2951 triads)
## parental triad expression estimates
PCmodel <- HEB[grepl("PC_",HEB$genotype),]
PCmodel <- PCmodel[PCmodel$group_id %in% gentoypecvpvals_pc$group_id,]

PC_Ppassfail <- PCmodel %>%
  dplyr::filter(genotype == "PC_P")
PC_Cpassfail <- PCmodel %>%
  dplyr::filter(genotype == "PC_C")

PC_P_Cpassfail <- PC_Ppassfail %>% 
  inner_join(PC_Cpassfail,by="group_id") %>%
  inner_join(gentoypecvpvals_pc,by='group_id')

PC_P_Cpassfail$gt_sig <- FALSE
PC_P_Cpassfail$gt_sig[PC_P_Cpassfail$adjustP < 0.0001] <- TRUE


# Now the same for PW

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

PWmodel <- HEB[grepl("PW_",HEB$genotype),]
PWmodel <- PWmodel[PWmodel$group_id %in% gentoypecvpvals_pw$group_id,]

PW_Ppassfail <- PWmodel %>%
  dplyr::filter(genotype == "PW_P")
PW_Wpassfail <- PWmodel %>%
  dplyr::filter(genotype == "PW_W")

PW_P_Wpassfail <- PW_Ppassfail %>% 
  inner_join(PW_Wpassfail,by="group_id") %>%
  inner_join(gentoypecvpvals_pw,by='group_id')

PW_P_Wpassfail$gt_sig <- FALSE
PW_P_Wpassfail$gt_sig[PW_P_Wpassfail$adjustP < 0.0001] <- TRUE
table(PW_P_Wpassfail$gt_sig)


# plot a barplot showing how many sig/nonsig overlap between the two F5 populations

PW_cv_lm_total_overlap <- unique(PW_P_Wpassfail[,c('group_id','gt_sig')])
PW_cv_lm_total_overlap$cross <- 'PW'
#Total PW triads ran through the CV~F5 line = 10,124

PC_cv_lm_total_overlap <- unique(PC_P_Cpassfail[,c('group_id','gt_sig')])
PC_cv_lm_total_overlap$cross <- 'PC'
#Total PC triads ran through the CV~F5 line = 11,109

PCPW_lm_total_overlap <- PC_cv_lm_total_overlap %>%
  inner_join(PW_cv_lm_total_overlap,by=c('group_id','gt_sig')) %>%
  dplyr::select(1:3)
colnames(PCPW_lm_total_overlap)[3] <- 'cross'
PCPW_lm_total_overlap$cross <- 'Overlap'


PCPW_lm_total_overlap_graph <- rbind(PW_cv_lm_total_overlap,PC_cv_lm_total_overlap,PCPW_lm_total_overlap)
PCPW_lm_total_overlap_graph$cross <- factor(PCPW_lm_total_overlap_graph$cross,levels = c('PC','PW','Overlap'))
PCPW_lm_total_overlap_graph$gt_sig[PCPW_lm_total_overlap_graph$gt_sig==TRUE] <- 'Significant effect\nof a genotype'
PCPW_lm_total_overlap_graph$gt_sig[PCPW_lm_total_overlap_graph$gt_sig==FALSE] <- 'No effect'

PCPW_cv_overlap <- PCPW_lm_total_overlap_graph
PCPW_cv_overlap$gt_sig <- gsub('Significant effect\nof a genotype','Significant effect of a genotype',PCPW_cv_overlap$gt_sig)
write.table(PCPW_cv_overlap,file = 'PCPW_lm_total_overlap_graph',sep = '\t',quote = F,row.names = F)

ggplot(PCPW_lm_total_overlap_graph,aes(cross)) + geom_bar_pattern(aes(pattern=gt_sig),color='black',fill='white',pattern_fill='black',
                                                                  pattern_spacing=0.03) +
  scale_pattern_manual(values = c('No effect' = 'none','Significant effect\nof a genotype'='stripe')) +
  scale_x_discrete(labels=c('PxC','PxW','Overlap')) +
  xlab('Cross') +
  ylab('Number of triads') +
  labs(pattern='') +
  theme_classic() +
  annotate("text",label = paste0(dim(PCPW_lm_total_overlap_graph[PCPW_lm_total_overlap_graph$gt_sig=='No effect' &
                                                                   PCPW_lm_total_overlap_graph$cross=='PC',])[1]),x=1,y=9000,size = 9) +
  annotate("text",label = paste0(dim(PCPW_lm_total_overlap_graph[PCPW_lm_total_overlap_graph$gt_sig=='No effect' &
                                                                   PCPW_lm_total_overlap_graph$cross=='PW',])[1]),x=2,y=8500,size = 9) +  
  annotate("text",label = paste0(dim(PCPW_lm_total_overlap_graph[PCPW_lm_total_overlap_graph$gt_sig=='No effect' &
                                                                   PCPW_lm_total_overlap_graph$cross=='Overlap',])[1]),x=3,y=5500,size = 9) +
  annotate(geom="rect",xmin=0.63,xmax=1.37,ymin=1600,ymax=2400,fill="white",color='black') +
  annotate(geom="rect",xmin=1.63,xmax=2.37,ymin=1400,ymax=2200,fill="white",color='black') +
  annotate(geom="rect",xmin=2.63,xmax=3.37,ymin=200,ymax=1000,fill="white",color='black') +
  annotate("text",label = paste0(dim(PCPW_lm_total_overlap_graph[PCPW_lm_total_overlap_graph$gt_sig=='Significant effect\nof a genotype' &
                                                                   PCPW_lm_total_overlap_graph$cross=='PC',])[1]),x=1,y=2000,size = 9) +
  annotate("text",label = paste0(dim(PCPW_lm_total_overlap_graph[PCPW_lm_total_overlap_graph$gt_sig=='Significant effect\nof a genotype' &
                                                                   PCPW_lm_total_overlap_graph$cross=='PW',])[1]),x=2,y=1800,size = 9) +
  annotate("text",label = paste0(dim(PCPW_lm_total_overlap_graph[PCPW_lm_total_overlap_graph$gt_sig=='Significant effect\nof a genotype' &
                                                                   PCPW_lm_total_overlap_graph$cross=='Overlap',])[1]),x=3,y=600,size = 9) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size=22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.position="top"
  )
ggsave('/Users/glombik/work/obj1_reanalysis/revised/cv_f6line_effect.png',device = 'png',height = 8,width = 10,dpi = 400,units = 'in',plot = last_plot())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/cv_f6line_effect.svg',device = 'svg',height = 8,width = 10,dpi = 400,units = 'in',plot = last_plot())

# check if chr position of genes has any pattern
homologies <- read.csv(file="/Users/glombik/work/obj1_reanalysis/homoeologs_1_1_1_synt_and_non_synt.csv")
lnghomologies <- homologies %>%
  pivot_longer(A:D,names_to = "genome")
colnames(lnghomologies)[8] <- c("gene.id")

PCPW_cv_overlap$group_id <- as.integer(PCPW_cv_overlap$group_id)
PCPW_cv_overlap_lng <- PCPW_cv_overlap %>%
  inner_join(lnghomologies,relationship = 'many-to-many')

gene_locations <- read.table('/Users/glombik/work/vcf_retry/geneIWGSC_v1.1_HC_20170706.gff3sorted')
colnames(gene_locations) <- c('chr','start','end','gene.id','strand')

PCPW_cv_overlap_lng <- PCPW_cv_overlap_lng %>%
  inner_join(gene_locations)

#Project it on a plot

ggplot(PCPW_cv_overlap_lng[PCPW_cv_overlap_lng$cross=='Overlap',],aes(x=start,fill=gt_sig)) + geom_histogram(alpha=0.8,bins=100) +
  facet_wrap(~chr,ncol=3) +
  theme_bw()
ggsave('/Users/glombik/work/obj1_reanalysis/revised/position_effect_cv_overlap_histogram.png',device = 'png',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

ggplot(PCPW_cv_overlap_lng[PCPW_cv_overlap_lng$gt_sig=='Significant effect of a genotype',],aes(x=start,fill=cross)) + 
  geom_histogram(alpha=0.8,bins=50) +
  facet_wrap(~chr,ncol=3) +
  theme_bw()
ggsave('/Users/glombik/work/obj1_reanalysis/revised/position_effect_cv_significant_triads_histogram.png',device = 'png',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

chr_len <- read.table('chr_len',header=T)
colnames(chr_len) <- c('chr','length')
chr_lenall <- chr_len
chr_lenall[2] <- 0
chr_lenall <- rbind(chr_lenall,chr_len)
data_with_lengths <- PCPW_cv_overlap_lng %>%
  inner_join(chr_len, by = "chr")

ggplot() +
  geom_blank(data=chr_len) +
  geom_density(data=data_with_lengths[data_with_lengths$gt_sig=='Significant effect of a genotype' &
                                        data_with_lengths$cross!='Overlap',],aes(x=start,fill=cross),alpha=0.5) +
  facet_wrap(~chr,ncol=3,scales = 'free') +
  scale_fill_manual(values = c('black','#009128ff')) +
  theme_bw()
ggsave('/Users/glombik/work/obj1_reanalysis/revised/position_effect_cv_significant_triads.png',device = 'png',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

ggplot() +
  geom_blank(data=chr_len) +
  geom_density(data=data_with_lengths[data_with_lengths$gt_sig=='Significant effect of a genotype' &
                                        data_with_lengths$cross!='Overlap',],aes(x=start,y=after_stat(count),fill=cross),alpha=0.5) +
  facet_wrap(~chr,ncol=3,scales = 'free') +
  scale_fill_manual(values = c('black','#009128ff')) +
  theme_bw()
ggsave('/Users/glombik/work/obj1_reanalysis/revised/position_effect_cv_significant_triads_rescaled_for_counts.png',device = 'png',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

sigdatawlengths <- data_with_lengths[data_with_lengths$gt_sig=='Significant effect of a genotype',]
table(sigdatawlengths$cross,sigdatawlengths$chr)

# Plot only chr A as the results are very similar across homoeologous chromosomes

chrAs_data_with_length <- data_with_lengths %>%
  filter(str_detect(chr,'A'))
chrAs_data_with_length$cross <- gsub('PC','PxC',chrAs_data_with_length$cross)
chrAs_data_with_length$cross <- gsub('PW','PxW',chrAs_data_with_length$cross)
chrA_len <- chr_len %>%
  filter(str_detect(chr,'A'))

ggplot() +
  geom_blank(data=chrA_len) +
  geom_density(data=chrAs_data_with_length[chrAs_data_with_length$gt_sig=='Significant effect of a genotype' &
                                             chrAs_data_with_length$cross!='Overlap',],aes(x=start,y=after_stat(count),fill=cross),alpha=0.5) +
  facet_wrap(~chr,ncol=1,scales = 'free') +
  scale_fill_manual(values = c("black", "#009128ff")) +
  theme_classic() +
  xlab('Position on chromosome') +
  ylab('Triad density') +
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 24),
    legend.text = element_text(size = 16),
    strip.text.x = element_text(size = 16),
    legend.position="top",
    legend.title = element_blank()
  )
ggsave('/Users/glombik/work/obj1_reanalysis/revised/onlyA_position_effect_cv_significant_triads_rescaled_for_counts.png',device = 'png',height = 12,width = 9,dpi = 400,units = 'in',plot = last_plot())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/onlyA_position_effect_cv_significant_triads_rescaled_for_counts.svg',device = 'svg',height = 12,width = 9,dpi = 400,units = 'in',plot = last_plot())



##### Now to combine it with bias distance results for a plot

#First I need distances between parents

## calculates distance matrix between parents
## this distance matrix is estimated by comparing first three columns to last three columns for each row separately
PC_all <- cbind(PC_P_Cpassfail[8:10],PC_P_Cpassfail[19:21])
PC_P_Cpassfail$distance <- apply(PC_all, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
PC_P_C_sum <- PC_P_Cpassfail[,c('group_id','gt_sig','distance')]
PC_P_C_sum$cross <-'PC'
write.table(PC_P_C_sum,'distPCpar.tsv',quote = F,row.names = F,sep = '\t')
# get parent triads with bias dist > 0.2 for check later
highdistPCpar <- PC_P_C_sum[PC_P_C_sum$distance>0.2,]
write.table(highdistPCpar,'highdistPCpar.tsv',quote = F,row.names = F,sep = '\t')
#how many % of triads display low bias dist? (<0.1)
length(unlist(PC_P_C_sum[PC_P_C_sum$distance<0.1,c('distance')]))/length(unlist(PC_P_C_sum[,c('distance')]))*100
#68.2 %

PW_all <- cbind(PW_P_Wpassfail[8:10],PW_P_Wpassfail[19:21])
PW_P_Wpassfail$distance <- apply(PW_all, 1, function(x) dist(matrix(x, nrow = 2, byrow = TRUE)))
PW_P_W_sum <- PW_P_Wpassfail[,c('group_id','gt_sig','distance')]
PW_P_W_sum$cross <- 'PW'
write.table(PW_P_W_sum,'distPWpar.tsv',quote = F,row.names = F,sep = '\t')

highdistPWpar <- PW_P_W_sum[PW_P_W_sum$distance>0.2,]
write.table(highdistPWpar,'highdistPWpar.tsv',quote = F,row.names = F,sep = '\t')
#how many % of triads display low bias dist? (<0.1)
length(unlist(PW_P_W_sum[PW_P_W_sum$distance<0.1,c('distance')]))/length(unlist(PW_P_W_sum[,c('distance')]))*100
#74.3 %

PCPW_sum <- rbind(PC_P_C_sum,PW_P_W_sum)
PCPW_sum$gt_sig[PCPW_sum$gt_sig==TRUE] <- 'Significant effect\nof a genotype'
PCPW_sum$gt_sig[PCPW_sum$gt_sig==FALSE] <- 'No effect'

PCPW_sum_corr_cross <- PCPW_sum
PCPW_sum_corr_cross$cross <- gsub('PC','PxC',PCPW_sum_corr_cross$cross)
PCPW_sum_corr_cross$cross <- gsub('PW','PxW',PCPW_sum_corr_cross$cross)

ggplot(PCPW_sum_corr_cross) + geom_density_pattern(aes(distance,pattern=gt_sig),color='black',fill='white',pattern_fill='black',
                                        pattern_spacing=0.03,alpha=0.3) +
  scale_pattern_manual(values = c('No effect' = 'none','Significant effect\nof a genotype'='stripe')) +
  facet_wrap(~cross,nrow = 2) +
  coord_cartesian(xlim = c(0,1.2)) +
  xlab('Distance between parents') +
  ylab('Scaled density') +
  theme_pubclean() +
  theme(axis.text = element_text(size = 22),
        axis.title = element_text(size=28),
        legend.title = element_blank(),
        legend.text = element_text(size = 28),
        strip.text.x = element_text(size = 22))
ggsave('/Users/glombik/work/obj1_reanalysis/revised/cv_f6line_distance_effect_zoomin.png',device = 'png',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
ggsave('/Users/glombik/work/obj1_reanalysis/revised/cv_f6line_distance_effect_zoomin.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())



#Perform a test on significance of bias distance between parents and significant effect of F5 line
head(PCPW_sum)
str(PCPW_sum)
PCPW_sum <- as.data.frame(PCPW_sum)
ks.test(PCPW_sum[PCPW_sum$cross=='PC' & PCPW_sum$gt_sig=='No effect',3],
        PCPW_sum[PCPW_sum$cross=='PC' & PCPW_sum$gt_sig=='Significant effect\nof a genotype',3])
ks.test(PCPW_sum[PCPW_sum$cross=='PW' & PCPW_sum$gt_sig=='No effect',3],
        PCPW_sum[PCPW_sum$cross=='PW' & PCPW_sum$gt_sig=='Significant effect\nof a genotype',3])

ks.test(distance~gt_sig,data = PCPW_sum[PCPW_sum$cross=='PW' & PCPW_sum$distance>0,])$p.value
ks.test(distance~gt_sig,data = PCPW_sum[PCPW_sum$cross=='PC' & PCPW_sum$distance>0,])$p.value


PCnoeffsample <- PCPW_sum[PCPW_sum$gt_sig=='No effect' & PCPW_sum$cross=='PC',] %>%
  sample_n(size = 2000)
PCsignifsample <- PCPW_sum[PCPW_sum$gt_sig=='Significant effect\nof a genotype' & PCPW_sum$cross=='PC',] %>%
  sample_n(size = 2000)
PCsample <- rbind(PCnoeffsample,PCsignifsample)
ks.test(distance~gt_sig,data = PCsample)$p.value

PWnoeffsample <- PCPW_sum[PCPW_sum$gt_sig=='No effect' & PCPW_sum$cross=='PW',] %>%
  sample_n(size = 2000)
PWsignifsample <- PCPW_sum[PCPW_sum$gt_sig=='Significant effect\nof a genotype' & PCPW_sum$cross=='PW',] %>%
  sample_n(size = 2000)
PWsample <- rbind(PWnoeffsample,PWsignifsample)
ks.test(distance~gt_sig,data = PWsample)
bothsample <- rbind(PCsample,PWsample)

# Plot the mean cv distribution in both populations
TPM <- read.table("/Users/glombik/work/obj1_reanalysis/TPMbeforeHEB_line193",header = T)
TPM$cross <- substr(TPM$genotype,1,2)
TPMcvgraph <- TPM %>%
  dplyr::select(c('group_id','sample','cv','cross')) %>%
  group_by(group_id,cross) %>%
  summarise(mean_cv_per_triad=mean(cv))


ggplot(TPMcvgraph) + geom_violin(aes(cross,mean_cv_per_triad,fill=cross),width=0.5,colour='white',alpha=0.2) + 
  geom_boxplot(aes(cross,mean_cv_per_triad,colour=cross),width=0.3,fill=alpha("white",alpha=0)) +
  # geom_jitter(aes(cross,mean_cv_per_triad,colour=cross),alpha=0.1,size=0.03,width = 0.1) +
  theme_classic() +
  xlab('F5 lines') +
  ylab('Mean CV per triad') +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 18),
    legend.text = element_text(size = 18),
    legend.title = element_text(size = 22)
  )
ggsave('/Users/glombik/work/obj1_reanalysis/revised/CV_variation_in_crosses.png',device = 'png',height = 8,width = 8,dpi = 400,units = 'in',plot = last_plot())

