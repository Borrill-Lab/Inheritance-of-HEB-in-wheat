# This script produces panels for figure 3 from intermediate data tables created during the analysis
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtern)
library(svglite)
library(ggpattern)
library(ggpubr)
library(rstatix)


# panel A
fig3_panelA <- read.table('/Users/glombik/work/Documents/my_articles/inheritance_2023/figure_scripts_github/fig3-panelA',header = T,sep = '\t')


#test for significant difference in log(SD) between associated and not associated homoeologs
wilcox_bothstdev <- fig3_panelA[fig3_panelA$logstdev!='-Inf',]
wilcox_bothstdev <- wilcox_bothstdev %>%
  group_by(cross) %>%
  wilcox_test(logstdev ~ classification) %>%
  add_significance() %>%
  add_xy_position(x='cross')
wilcox_bothstdev

ggplot(fig3_panelA,aes(cross,logstdev)) +
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
ggsave(plot = last_plot(),filename = 'fig3_panelA.svg',device = 'svg',dpi = 400,width = 10,height = 10,units = "in")

# panel B
fig3_panelB <- read.table('/Users/glombik/work/Documents/my_articles/inheritance_2023/figure_scripts_github/fig3-panelB',header = T,sep = '\t')

ggplot(fig3_panelB,aes(`To.C`,`To.P`)) + 
  geom_point(data=fig3_panelB,size=5,aes(group=inherited_genotype_collapsed,fill=inherited_genotype_collapsed),pch=21,alpha=0.8,color='black',show.legend = F) +
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

ggsave('fig3_panelB.svg',device = 'svg',height = 8,
       width = 8,dpi = 400,units = 'in',plot = last_plot())

# panel C - not final - finishing touches done in Inkscape
fig3_panelC <- read.table('/Users/glombik/work/Documents/my_articles/inheritance_2023/figure_scripts_github/fig3-panelC',header = T,sep = '\t')
fig3parentgroup = c("Charger = P2","Paragon = P1")

ggtern(fig3_panelC,aes(A_tpm,D_tpm,B_tpm),group=inherited_genotype_collapsed) +
  geom_Tline(Tintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_Rline(Rintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_Lline(Lintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_mask() +
  geom_point(data = fig3_panelC %>% filter(!inherited_genotype_collapsed %in% fig3parentgroup),aes(fill=inherited_genotype_collapsed),pch=21,colour="black",size=4,alpha=0.8) +
  geom_point(data = fig3_panelC %>% filter(inherited_genotype_collapsed %in% fig3parentgroup),aes(fill=inherited_genotype_collapsed),pch=21,colour="black",size=8,alpha=0.8) +
  geom_point(data = fig3_panelC %>% filter(!inherited_genotype_collapsed %in% fig3parentgroup),aes(shape=associated_homoeologue),alpha=0) +
  scale_fill_manual(values = c('#fff582ff','#DC4462','#7F9596','#efa9a9ff')) +
  labs(x="A",y="D",z="B") +
  theme_gray() +
  theme(axis.title = element_text(size = 15),axis.text = element_text(size = 10))

ggsave(filename = paste0("fig3_panelC.svg"),device = "svg",dpi = 400,
       plot = last_plot(),width = 8,height = 8,units = "in")
