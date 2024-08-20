# This script produces panels for figure 1 from intermediate data tables created during the analysis
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtern)
library(svglite)
library(ggpattern)
library(ggpubr)
# Panel A - base with an example triad expression based on CV difference between replicates - finishing touches done in Inkscape
fig1_panelA <- read.table('/Users/glombik/work/Documents/my_articles/inheritance_2023/figure_scripts_github/fig1-panelA',header = T)

ggtern(fig1_panelA,aes(A_tpm,D_tpm,B_tpm)) +
  geom_Tline(Tintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_Rline(Rintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_Lline(Lintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_mask() +
  geom_point(colour="black",size=3,alpha=0.8,fill='white',pch=21) +
  labs(x="A",y="D",z="B") +
  theme_gray() +
  theme(axis.title = element_text(size = 15),axis.text = element_text(size = 10))
ggsave('fig1_panelA.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

# Panel B

fig1_panelB <- read.table('/Users/glombik/work/Documents/my_articles/inheritance_2023/figure_scripts_github/fig1-panelB',header = T)
fig1_panelB$Categories <- as.character(fig1_panelB$Categories)

ggplot(fig1_panelB,aes(perc,countperc,color=cross,linetype=Categories)) +
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

ggsave('fig1-panelB.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

# Panel C
fig1_panelC <- read.table('/Users/glombik/work/Documents/my_articles/inheritance_2023/figure_scripts_github/fig1-panelC',header = T,sep = '\t')
fig1_panelC$gt_sig <- gsub('Significant effect of a genotype','Significant effect\nof a genotype',fig1_panelC$gt_sig)

ggplot(fig1_panelC,aes(cross)) + geom_bar_pattern(aes(pattern=gt_sig),color='black',fill='white',pattern_fill='black',
                                                                  pattern_spacing=0.03) +
  scale_pattern_manual(values = c('No effect' = 'none','Significant effect\nof a genotype'='stripe')) +
  scale_x_discrete(labels=c('PxC','PxW','Overlap')) +
  xlab('Cross') +
  ylab('Number of triads') +
  labs(pattern='') +
  theme_classic() +
  annotate("text",label = paste0(dim(fig1_panelC[fig1_panelC$gt_sig=='No effect' &
                                                   fig1_panelC$cross=='PC',])[1]),x=1,y=9000,size = 9) +
  annotate("text",label = paste0(dim(fig1_panelC[fig1_panelC$gt_sig=='No effect' &
                                                   fig1_panelC$cross=='PW',])[1]),x=2,y=8500,size = 9) +  
  annotate("text",label = paste0(dim(fig1_panelC[fig1_panelC$gt_sig=='No effect' &
                                                   fig1_panelC$cross=='Overlap',])[1]),x=3,y=5500,size = 9) +
  annotate(geom="rect",xmin=0.63,xmax=1.37,ymin=1600,ymax=2400,fill="white",color='black') +
  annotate(geom="rect",xmin=1.63,xmax=2.37,ymin=1400,ymax=2200,fill="white",color='black') +
  annotate(geom="rect",xmin=2.63,xmax=3.37,ymin=200,ymax=1000,fill="white",color='black') +
  annotate("text",label = paste0(dim(fig1_panelC[fig1_panelC$gt_sig=='Significant effect\nof a genotype' &
                                                   fig1_panelC$cross=='PC',])[1]),x=1,y=2000,size = 9) +
  annotate("text",label = paste0(dim(fig1_panelC[fig1_panelC$gt_sig=='Significant effect\nof a genotype' &
                                                   fig1_panelC$cross=='PW',])[1]),x=2,y=1800,size = 9) +
  annotate("text",label = paste0(dim(fig1_panelC[fig1_panelC$gt_sig=='Significant effect\nof a genotype' &
                                                   fig1_panelC$cross=='Overlap',])[1]),x=3,y=600,size = 9) +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size=22),
        legend.title = element_text(size = 22),
        legend.text = element_text(size = 18),
        legend.position="top"
  )
ggsave('fig1-panelC.svg',device = 'svg',height = 8,width = 10,dpi = 400,units = 'in',plot = last_plot())

# Panel D - Plot only chr A as the results are very similar across homoeologous chromosomes
fig1_panelD <- read.table('/Users/glombik/work/Documents/my_articles/inheritance_2023/figure_scripts_github/fig1-panelD',header = T,sep = '\t')

chr_len <- read.table('chr_len',header=T)
colnames(chr_len) <- c('chr','length')
chr_lenall <- chr_len
chr_lenall[2] <- 0
chr_lenall <- rbind(chr_lenall,chr_len)
fig1_data_with_lengths <- fig1_panelD %>%
  inner_join(chr_len, by = "chr")
chrA_len <- chr_len %>%
  filter(str_detect(chr,'A'))

fig1_chrAs_data_with_length <- fig1_data_with_lengths %>%
  filter(str_detect(chr,'A'))
fig1_chrAs_data_with_length$cross <- gsub('PC','PxC',fig1_chrAs_data_with_length$cross)
fig1_chrAs_data_with_length$cross <- gsub('PW','PxW',fig1_chrAs_data_with_length$cross)

ggplot() +
  geom_blank(data=chrA_len) +
  geom_density(data=fig1_chrAs_data_with_length[fig1_chrAs_data_with_length$gt_sig=='Significant effect of a genotype' &
                                                  fig1_chrAs_data_with_length$cross!='Overlap',],aes(x=start,y=after_stat(count),fill=cross),alpha=0.5) +
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
ggsave('fig1-panelD.svg',device = 'svg',height = 12,width = 9,dpi = 400,units = 'in',plot = last_plot())

# Panel E
fig1_panelE <- read.table('/Users/glombik/work/Documents/my_articles/inheritance_2023/figure_scripts_github/fig1-panelE',header = T,sep = '\t')
fig1_panelE$gt_sig <- gsub('Significant effect of a genotype','Significant effect\nof a genotype',fig1_panelE$gt_sig)

ggplot(fig1_panelE) + geom_density_pattern(aes(distance,pattern=gt_sig),color='black',fill='white',pattern_fill='black',
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
ggsave('fig1_panelE.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
