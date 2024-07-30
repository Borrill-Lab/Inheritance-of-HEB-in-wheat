# Marek Glombik
# This script does statistics about bias distance between parental lines and F5 lines and produces panels for figures 1-2

library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggtern)
library(LSD)
library(ggExtra)
library(gridExtra)
library(gplots)

F6_PCwhole <- read.table('/Users/glombik/work/obj1_reanalysis/F6_PCwhole',header = T,sep = '\t')
PCdistcheck <- F6_PCwhole %>%
  pivot_wider(values_from = c('distance_to_parent','cv_parent','tpm_parent'),names_from = which_parent)

heatscatter(PCdistcheck$`distance_to_parent_To C`,PCdistcheck$`distance_to_parent_To P`,xlab = "distance to paragon in P x C",ylab="distance to charger in P x C")
abline(0,1)

heathist(PCdistcheck$`distance_to_parent_To C`)
heathist(PCdistcheck$`distance_to_parent_To P`)

mean(PCdistcheck$`distance_to_parent_To C`)
mean(PCdistcheck$`distance_to_parent_To P`)

F6_PWwhole <- read.table('/Users/glombik/work/obj1_reanalysis/F6_PWwhole',header = T,sep = '\t')
PWdistcheck <- F6_PWwhole %>%
  pivot_wider(values_from = c('distance_to_parent','cv_parent','tpm_parent'),names_from = which_parent)

heatscatter(PWdistcheck$`distance_to_parent_To W`,PWdistcheck$`distance_to_parent_To P`)
heathist(PWdistcheck$`distance_to_parent_To W`)
heathist(PWdistcheck$`distance_to_parent_To P`)

mean(PWdistcheck$`distance_to_parent_To W`)
mean(PWdistcheck$`distance_to_parent_To P`)

# Mean distance in both populations is 0.11 and median 0.07...

PChist <- F6_PCwhole[,c('distance_to_parent','which_parent')]
PChist$which_parent <- gsub('$',' PxC',PChist$which_parent)
PWhist <- F6_PWwhole[,c('distance_to_parent','which_parent')]
PWhist$which_parent <- gsub('$',' PxW',PWhist$which_parent)

PCPWhist <- rbind(PChist,PWhist)
library(ggplot2)

mean_data <- PCPWhist %>%
  group_by(which_parent) %>%
  summarise(mean_distance = mean(distance_to_parent))


ggplot(PCPWhist,aes(x=distance_to_parent,fill=which_parent,color=which_parent)) +
  geom_density(position = 'stack',show.legend = F) +
  geom_vline(aes(xintercept=mean(distance_to_parent)),color='black',linetype='dashed',show.legend = F,linewidth=2) +
  geom_text(data = mean_data, aes(x = 0.4, y = 7, color = which_parent, label = paste0('mean = ', round(mean_distance, digits = 1))),
            position = 'stack',size=10,show.legend = F) +
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  scale_color_manual(values = c("#999999", "#E69F00", "#56B4E9", "#009E73")) +
  ylab(label = 'number of F5 lines') +
  xlab(label = 'bias distance to parent') +
  theme_classic() +
  theme(axis.title = element_text(size=23),
        axis.text = element_text(size=16)) +
  annotate("label", x = 0.9, y = 28, label = "Bias distance to Charger (PxC)", label.size = 2.0,size=10, colour="white",fill='#999999') +
  annotate("label", x = 0.9, y = 21, label = "Bias distance to Paragon (PxC)", label.size = 2.0, cex=10, fill="#E69F00",colour='white') +
  annotate("label", x = 0.9, y = 14, label = "Bias distance to Watkins (PxW)", label.size = 2.0, cex=10, fill="#56B4E9",colour='white') +
  annotate("label", x = 0.9, y = 7, label = "Bias distance to Paragon (PxW)", label.size = 2.0, cex=10, fill="#009E73",colour='white')
  
  
ggsave(plot = last_plot(),filename = 'bias_dist_hist.png',device = 'png',dpi = 400,width = 14,height = 10,units = "in")

PChistbelow02 <- PChist[PChist$distance_to_parent<0.12,]
PWhistbelow02 <- PWhist[PWhist$distance_to_parent<0.12,]

# In triads that are biased towards one or the other parent, we want to check the proportion of F5 lines towards each parent

PCbiased <- read.table('/Users/glombik/work/obj1_reanalysis/PC_P_C_sum_biased',header = T,sep = '\t')
table(PCbiased$divergence_in_F6)


cntPCdist <- PCdistcheck %>%
  dplyr::group_by(group_id) %>%
  dplyr::summarize(
    num_genotypes_followP_PC = sum(`distance_to_parent_To C` > 0.2 & `distance_to_parent_To P` < 0.1, na.rm = TRUE),
    num_genotypes_followC_PC = sum(`distance_to_parent_To P` > 0.2 & `distance_to_parent_To C` < 0.1, na.rm = TRUE)
  )

# now take those which have > 15 in both

cntPCdistbothover15 <- cntPCdist[cntPCdist$num_genotypes_followP_PC>15 & cntPCdist$num_genotypes_followC_PC>15,]
cntPCdistbothover15$rowgroup <- rownames(cntPCdistbothover15)

cntPCdistbothover15 <- cntPCdistbothover15 %>%
  pivot_longer(!c(group_id,rowgroup),names_to = 'parent',values_to = 'genotype_count' )

ggplot(cntPCdistbothover15) +
  geom_boxplot(aes(x=parent,y=genotype_count,fill=parent)) +
  geom_point(aes(x=parent,y=genotype_count),position = 'jitter',size=1)

mean_P <- mean(subset(cntPCdistbothover15, parent == 'num_genotypes_followP_PC')$genotype_count)
mean_C <- mean(subset(cntPCdistbothover15, parent == 'num_genotypes_followC_PC')$genotype_count)


PWbiased <- read.table('/Users/glombik/work/obj1_reanalysis/PW_P_W_sum_biased',header = T,sep = '\t')
table(PWbiased$divergence_in_F6)

cntPWdist <- PWdistcheck %>%
  dplyr::group_by(group_id) %>%
  dplyr::summarize(
    num_genotypes_followP_PW = sum(`distance_to_parent_To W` > 0.2 & `distance_to_parent_To P` < 0.1, na.rm = TRUE),
    num_genotypes_followW_PW = sum(`distance_to_parent_To P` > 0.2 & `distance_to_parent_To W` < 0.1, na.rm = TRUE)
  )

# now take those which have > 15 in both

cntPWdistbothover15 <- cntPWdist[cntPWdist$num_genotypes_followP_PW>15 & cntPWdist$num_genotypes_followW_PW>15,]
cntPWdistbothover15$rowgroup <- rownames(cntPWdistbothover15)

cntPWdistbothover15 <- cntPWdistbothover15 %>%
  pivot_longer(!c(group_id,rowgroup),names_to = 'parent',values_to = 'genotype_count' )

ggplot(cntPWdistbothover15) +
  geom_boxplot(aes(x=parent,y=genotype_count,fill=parent)) +
  geom_point(aes(x=parent,y=genotype_count),position = 'jitter',size=1)

mean_P <- mean(subset(cntPWdistbothover15, parent == 'num_genotypes_followP_PW')$genotype_count)
mean_W <- mean(subset(cntPWdistbothover15, parent == 'num_genotypes_followW_PW')$genotype_count)

# Put them together

cntdist_both_crosses <- rbind(cntPCdistbothover15,cntPWdistbothover15)
cntdist_both_crosses$parent <- factor(cntdist_both_crosses$parent,levels = c('num_genotypes_followP_PC','num_genotypes_followC_PC','num_genotypes_followP_PW',
                                                 'num_genotypes_followW_PW'))
ggplot(cntdist_both_crosses) +
  geom_boxplot(aes(x=parent,y=genotype_count,fill=parent),show.legend=F) +
  geom_point(aes(x=parent,y=genotype_count),position = 'jitter',size=1,alpha=0.4) +
  scale_fill_manual(values = c("#999999", "#E69F00", "#009E73", "#56B4E9")) +
  ylab(expression(label = 'Number of F'[5]~' lines following\n parental HEB in a triad')) +
  xlab(label = 'Parent') +
  ylim(0,35) +
  theme_classic() +
  theme(axis.title = element_text(size=38),
        axis.text = element_text(size=32)) +
  scale_x_discrete(labels=c('num_genotypes_followP_PC'='Paragon\n(PxC)','num_genotypes_followC_PC'='Charger\n(PxC)',
                            'num_genotypes_followW_PW'='Watkins\n(PxW)','num_genotypes_followP_PW'='Paragon\n(PxW)'))

ggsave(plot = last_plot(),filename = 'bias_dist_boxplot.svg',device = 'svg',dpi = 400,width = 14,height = 10,units = "in")


# Make the heatscatter plots better based on Philippa's suggestions 
# First try to plot it with ggplot
heatscatter(PCdistcheck$`distance_to_parent_To C`,PCdistcheck$`distance_to_parent_To P`,xlab = "distance to paragon in P x C",ylab="distance to charger in P x C")

ggplot(PCdistcheck,aes(`distance_to_parent_To C`,`distance_to_parent_To P`)) + 
  # geom_point(fill='grey',color='grey') +
  geom_bin2d(bins = 200, aes(fill = stat(density)),show.legend = F) +
  scale_fill_gradientn(colors = c('grey','blue','red','yellow'),
                       breaks=c(0.004,0.008,0.010,0.012),
                       labels=c('low','mid','midhigh','high')) +
  xlim(0,1.2) +
  ylim(0,1.2) +
  theme_classic()


F6_PWwhole <- read.table('/Users/glombik/work/obj1_reanalysis/F6_PWwhole',header = T,sep = '\t')
PWdistcheck <- F6_PWwhole %>%
  pivot_wider(values_from = c('distance_to_parent','cv_parent','tpm_parent'),names_from = which_parent)

heatscatter(PWdistcheck$`distance_to_parent_To W`,PWdistcheck$`distance_to_parent_To P`)

ggplot(PWdistcheck,aes(`distance_to_parent_To W`,`distance_to_parent_To P`)) + 
  # geom_point(fill='grey',color='grey') +
  geom_bin2d(bins = 200, aes(fill = stat(density)),show.legend = F) +
  scale_fill_gradientn(colors = c('grey','blue','red','yellow'),
                       breaks=c(0.004,0.008,0.010,0.012),
                       labels=c('low','mid','midhigh','high')) +
  xlim(0,1.2) +
  ylim(0,1.2) +
  theme_classic()

# Plots can be made, now colour parts of the graph by category of bias dist with annotate rectangles
# category 1 is dist < 0.1 to C and dist > 0.2 to P
# category 2 is dist < 0.1 to P and dist > 0.2 to C
# category 3 is dist > 0.2 to C and dist > 0.2 to P
# category 4 is everything else


ggplot(PCdistcheck,aes(`distance_to_parent_To C`,`distance_to_parent_To P`)) + 
  # geom_point(fill='grey',color='grey') +
  geom_bin2d(bins = 200, aes(fill = stat(density)),show.legend = F,) +
  scale_fill_gradientn(colors = c('grey','blue','red','yellow'),
                       breaks=c(0.004,0.008,0.010,0.012),
                       labels=c('low','mid','midhigh','high')) +
  annotate('rect',xmin=0,xmax=0.1,ymin=0.2,ymax=1.2,fill='yellow',alpha=0.3) +
  annotate('rect',xmin=0.2,xmax=1.2,ymin=0,ymax=0.1,fill='yellow',alpha=0.3) +
  annotate('rect',xmin=0.2,xmax=1.2,ymin=0.2,ymax=1.2,fill='red',alpha=0.2) +
  annotate('text',x=0.05,y=0.8,label='1',size=9) +
  annotate('text',x=0.8,y=0.05,label='2',size=9) +
  annotate('text',x=0.8,y=0.8,label='3',size=9) +
  annotate('text',x=0.8,y=0.8,label='4',size=9) +
  xlim(0,1.2) +
  ylim(0,1.2) +
  scale_x_continuous(breaks = seq(0,1.2,by=0.2)) +
  scale_y_continuous(breaks = seq(0,1.2,by=0.2)) +
  theme_classic()


# Now I need to add a column category to separate them otherwise it will just plot all points goruped together
library(ggExtra)


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




ggplot(PWdistcheck,aes(`distance_to_parent_To W`,`distance_to_parent_To P`)) + 
  geom_point(data=PWdistcheck,size=0.3,aes(group=category,color=category),show.legend = F,alpha=0) +
  scale_color_manual(values = c('#999999','#56B4E9','#009E73','#F0E442','#E69F00')) +
  geom_bin2d(bins = 200, aes(fill = stat(density)),show.legend = F,) +
  scale_fill_gradientn(colors = c('grey','blue','red','yellow'),
                       breaks=c(0.004,0.008,0.010,0.012),
                       labels=c('low','mid','midhigh','high')) +
  annotate('rect',xmin=0,xmax=0.1,ymin=0.2,ymax=1.2,fill='#56B4E9',alpha=0.15,colour='black', linetype='dashed') +
  annotate('rect',xmin=0.2,xmax=1.2,ymin=0,ymax=0.1,fill='#009E73',alpha=0.15,colour='black', linetype='dashed') +
  annotate('rect',xmin=0.2,xmax=1.2,ymin=0.2,ymax=1.2,fill='#F0E442',alpha=0.15,colour='black', linetype='dashed') +
  annotate('rect',xmin=0,xmax=0.2,ymin=0,ymax=0.2,fill='#E69F00',alpha=0.15,colour='black', linetype='dashed') +
  annotate('text',x=0.05,y=0.8,label='1',size=13) +
  annotate('text',x=0.8,y=0.05,label='2',size=13) +
  annotate('text',x=0.8,y=0.8,label='3',size=13) +
  annotate('text',x=0.1,y=0.1,label='4',size=13) +
  xlab('Bias distance to Watkins parent') +
  ylab('Bias distance to Paragon parent') +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size=24),
        axis.text = element_text(size = 20))
ggsave(filename = paste0("/Users/glombik/work/Documents/my_articles/inheritance_2023/bias_dist_PW_adj_no_margins",".png"),device = "png",dpi = 400,
       plot = last_plot(),width = 10,height = 8,units = "in")
ggsave(filename = paste0("/Users/glombik/work/Documents/my_articles/inheritance_2023/bias_dist_PW_adj_no_margins",".svg"),device = "svg",dpi = 400,
       plot = last_plot(),width = 10,height = 8,units = "in")

ggplot(PCdistcheck,aes(`distance_to_parent_To C`,`distance_to_parent_To P`)) + 
  geom_point(data=PCdistcheck,size=0.3,aes(group=category,color=category),show.legend = F,alpha=0) +
  scale_color_manual(values = c('#999999','#56B4E9','#009E73','#F0E442','#E69F00')) +
  geom_bin2d(bins = 200, aes(fill = stat(density)),show.legend = F,) +
  scale_fill_gradientn(colors = c('grey','blue','red','yellow'),
                       breaks=c(0.004,0.008,0.010,0.012),
                       labels=c('low','mid','midhigh','high')) +
  annotate('rect',xmin=0,xmax=0.1,ymin=0.2,ymax=1.2,fill='#56B4E9',alpha=0.15,colour='black', linetype='dashed') +
  annotate('rect',xmin=0.2,xmax=1.2,ymin=0,ymax=0.1,fill='#009E73',alpha=0.15,colour='black', linetype='dashed') +
  annotate('rect',xmin=0.2,xmax=1.2,ymin=0.2,ymax=1.2,fill='#F0E442',alpha=0.15,colour='black', linetype='dashed') +
  annotate('rect',xmin=0,xmax=0.2,ymin=0,ymax=0.2,fill='#E69F00',alpha=0.15,colour='black', linetype='dashed') +
  annotate('text',x=0.05,y=0.8,label='1',size=13) +
  annotate('text',x=0.8,y=0.05,label='2',size=13) +
  annotate('text',x=0.8,y=0.8,label='3',size=13) +
  annotate('text',x=0.1,y=0.1,label='4',size=13) +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) +
  xlab('Bias distance to Charger parent') +
  ylab('Bias distance to Paragon parent') +
  theme_classic() +
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size=24),
        axis.text = element_text(size = 20))

ggsave(filename = paste0("/Users/glombik/work/Documents/my_articles/inheritance_2023/bias_dist_PC_adj_no_margins",".png"),device = "png",dpi = 400,
       plot = last_plot(),width = 10,height = 8,units = "in")
ggsave(filename = paste0("/Users/glombik/work/Documents/my_articles/inheritance_2023/bias_dist_PC_adj_no_margins",".svg"),device = "svg",dpi = 400,
       plot = last_plot(),width = 10,height = 8,units = "in")


pctriads <- unique(PCdistcheck[,c('group_id')])
pwtriads <- unique(PWdistcheck[,c('group_id')])
jointriads <- pctriads %>%
  inner_join(pwtriads)


# Make a summary bar graph with the categories 1-4
PCsumbar <- PCdistcheck[,c('group_id','category')]
PWsumbar <- PWdistcheck[,c('group_id','category')]

PCsubmartab <- as.data.frame(table(PCsumbar$group_id,PCsumbar$category))
PWsubmartab <- as.data.frame(table(PWsumbar$group_id,PWsumbar$category))
colnames(PCsubmartab) <- c('group_id','category','counts')
colnames(PWsubmartab) <- c('group_id','category','counts')
PCsubmartab$cat_over_15 <- NA
PWsubmartab$cat_over_15 <- NA
PCsubmartab$cat_over_15 <- ifelse(PCsubmartab$counts > 15, paste(PCsubmartab$category), PCsubmartab$cat_over_15)
PWsubmartab$cat_over_15 <- ifelse(PWsubmartab$counts > 15, paste(PWsubmartab$category), PWsubmartab$cat_over_15)
PCgraph_sumbartab <- PCsubmartab[!is.na(PCsubmartab$cat_over_15),]
PWgraph_sumbartab <- PWsubmartab[!is.na(PWsubmartab$cat_over_15),]

pcna_graph <- PCsubmartab[is.na(PCsubmartab$cat_over_15),]

PCgraph_sumbartab$cross <- 'PxC'
PWgraph_sumbartab$cross <- 'PxW'
overlap_sumbartab <- PCgraph_sumbartab %>%
  inner_join(PWgraph_sumbartab,by=c('group_id','cat_over_15'))
overlap_sumbartab <- overlap_sumbartab[,1:5]
colnames(overlap_sumbartab) <- c('group_id','category','counts','cat_over_15','cross')
overlap_sumbartab$cross <- 'Overlap'

#how many DFO/DFB triads overlap?
length(unlist(unique(overlap_sumbartab[overlap_sumbartab$cat_over_15!='Conserved' &
                                         overlap_sumbartab$cat_over_15!='Uncategorised',c('group_id')])))


bothgraph_sumbartab <- rbind(PCgraph_sumbartab,PWgraph_sumbartab,overlap_sumbartab)
bothgraph_sumbartab$cat_over_15 <- factor(bothgraph_sumbartab$cat_over_15,levels=c('Uncategorised','DFO_b','DFO_a','DFB','Conserved'))
bothgraph_sumbartab$cross <- factor(bothgraph_sumbartab$cross,levels=c('PxC','PxW','Overlap'))


ggplot(bothgraph_sumbartab) + geom_bar(aes(cat_over_15,fill=cross),colour='black',position = position_dodge2(reverse=T)) +
  xlab('Category') +
  ylab('Triad count') +
  scale_fill_manual(values = c('black','#009128ff','yellow')) +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 30),
        axis.text.x = element_text(angle=45),
        axis.title = element_text(size = 34),
        axis.title.y = element_text(angle = 90),
        legend.text = element_text(size = 30),
        legend.title = element_blank())
ggsave(plot = last_plot(),filename = 'bias_dist_categories_flip.png',device = 'png',dpi = 400,width = 10,height = 10,units = "in")
ggsave(plot = last_plot(),filename = 'bias_dist_categories_flip.svg',device = 'svg',dpi = 400,width = 10,height = 10,units = "in")

PCPWgraph_sumbartab <- rbind(PCgraph_sumbartab,PWgraph_sumbartab)
write.table(PCPWgraph_sumbartab,'PCPWbiasdist_sum_categories.tsv',row.names = F,quote = F,sep = '\t')


### Now do a heatmap of how many % of F5 lines per triad are in eaach categ
heatPCsumbartab <- PCsubmartab[,1:3]
heatPCsumbartab$perc <- (heatPCsumbartab$counts/50)*100

hheatPC <- heatPCsumbartab %>%
  dplyr::select(group_id,category,perc) %>%
  pivot_wider(names_from = category,values_from = perc)
hheatPC <- as.data.frame(hheatPC)
rownames(hheatPC) <- hheatPC$group_id
hheatPC <- hheatPC[,-1]

transposed_data <- t(hheatPC)
pdf('heat_triads_distance_PC.pdf',width = 12,height = 9)
svg('heat_triads_distance_PC.svg',width = 12,height = 9)
png('heat_triads_distance_PC.png',width = 12,height = 9,units = 'in',res=1000)
heatmap.2(
  transposed_data,
  col = colorRampPalette(c("#F1F0EF" ,"black"))(100), # Color palette
  main = "",
  trace = 'none',
  cexRow = 0.5,
  Rowv = F,
  labCol = FALSE,
  srtRow = 45,
  dendrogram = 'none',
  key.xlab = 'Percentage',
  key.title = '',
  sepwidth = c(0,0),
  density.info = 'none'
)
dev.off()


heatPWsumbartab <- PWsubmartab[,1:3]
heatPWsumbartab$perc <- (heatPWsumbartab$counts/50)*100

hheatPW <- heatPWsumbartab %>%
  dplyr::select(group_id,category,perc) %>%
  pivot_wider(names_from = category,values_from = perc)
hheatPW <- as.data.frame(hheatPW)
rownames(hheatPW) <- hheatPW$group_id
hheatPW <- hheatPW[,-1]

transposed_dataPW <- t(hheatPW)
pdf('heat_triads_distance_PW.pdf',width = 12,height = 9)
svg('heat_triads_distance_PW.svg',width = 12,height = 9)
png('heat_triads_distance_PW.png',width = 12,height = 9,units = 'in',res=1000)
heatmap.2(
  transposed_dataPW,
  col = colorRampPalette(c("#F1F0EF" ,"#009128ff"))(100), # Color palette
  main = "",
  trace = 'none',
  cexRow = 0.5,
  Rowv = F,
  labCol = FALSE,
  srtRow = 45,
  dendrogram = 'none',
  key.xlab = 'Percentage',
  key.ylab = '',
  key.title = '',
  sepwidth = c(0,0),
  density.info = 'none'
  )
dev.off()



# Plot the distribution of bias distance in parental lines in both populations for DFB category
PCPWgraph_sumbartab <- read.table('/Users/glombik/work/vcf_retry/PCPWbiasdist_sum_categories.tsv',header = T)
head(PCPWgraph_sumbartab)
head(distPCpar)

distPCpar <- read.table('/Users/glombik/work/obj1_reanalysis/distPCpar.tsv',header = T)
distPWpar <- read.table('/Users/glombik/work/obj1_reanalysis/distPWpar.tsv',header = T)

dfbPC <- PCPWgraph_sumbartab[PCPWgraph_sumbartab$cross=='PxC' & PCPWgraph_sumbartab$cat_over_15=='DFB',]
dfbPW <- PCPWgraph_sumbartab[PCPWgraph_sumbartab$cross=='PxW' & PCPWgraph_sumbartab$cat_over_15=='DFB',]

dfbPC <- dfbPC %>%
  inner_join(distPCpar,by='group_id')
dfbPW <- dfbPW %>%
  inner_join(distPWpar,by='group_id')
dfbboth <- rbind(dfbPC,dfbPW)

ggplot(dfbboth,aes(cross.x,distance,color=cross.x)) +
  geom_violin(show.legend = F) +
  geom_point(position = 'jitter',size=1,alpha=0.6,show.legend = F) +
  scale_color_manual(values = c('black','#009128ff')) +
  theme_classic() +
  xlab('Cross') +
  ylab('Bias distance between\ntriads in parental lines') +
  theme(axis.text = element_text(size = 32),
        axis.title = element_text(size=38))
ggsave('/Users/glombik/work/obj1_reanalysis/revised/DFB_bias_dist_parents.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())


length(unlist(dfbPC[dfbPC$distance>0.2,c('distance')]))
#49 out of 56 = 87.5 %
length(unlist(dfbPW[dfbPW$distance>0.2,c('distance')]))
#50 out of 309 = 16.2 %
length(unlist(dfbPW[dfbPW$distance<0.2,c('distance')]))
#259 out of 309 = 83.8 %

