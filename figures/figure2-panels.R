# This script produces panels for figure 2 from intermediate data tables created during the analysis
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggtern)
library(svglite)
library(ggpattern)
library(ggpubr)
library(gplots)

# panel A
fig2_panelA <- read.table('fig2-panelA',header = T,sep = '\t')


ggplot(fig2_panelA,aes(`distance_to_parent_To.C`,`distance_to_parent_To.P`)) + 
  geom_point(data=fig2_panelA,size=0.3,aes(group=category,color=category),show.legend = F,alpha=0) +
  scale_color_manual(values = c('#999999','#56B4E9','#009E73','#F0E442','#E69F00')) +
  geom_bin2d(bins = 200, aes(fill = stat(density)),show.legend = F,) +
  scale_fill_gradientn(colors = c('grey','blue','red','yellow'),
                       breaks=c(0.004,0.008,0.010,0.012),
                       labels=c('low','mid','midhigh','high')) +
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
ggsave('fig2_panelA.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

# panel B
fig2_panelB <- read.table('fig2-panelB',header = T,sep = '\t')

ggplot(fig2_panelB,aes(`distance_to_parent_To.W`,`distance_to_parent_To.P`)) + 
  geom_point(data=fig2_panelB,size=0.3,aes(group=category,color=category),show.legend = F,alpha=0) +
  scale_color_manual(values = c('#999999','#56B4E9','#009E73','#F0E442','#E69F00')) +
  geom_bin2d(bins = 200, aes(fill = stat(density)),show.legend = F,) +
  scale_fill_gradientn(colors = c('grey','blue','red','yellow'),
                       breaks=c(0.004,0.008,0.010,0.012),
                       labels=c('low','mid','midhigh','high')) +
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
  xlab('Bias distance to Watkins parent') +
  ylab('Bias distance to Paragon parent') +
  scale_x_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) +
  scale_y_continuous(breaks = seq(0, 1.2, by = 0.2), limits = c(0, 1.2)) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 24),
        axis.title.y = element_text(size=24),
        axis.text = element_text(size = 25))
ggsave('fig2_panelB.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())


# panel C - for this panel only empty triangular plot was used - finishing touches were done in Inkscape
emptycv <- data.frame(
  A_tpm = c(2, 5, 3),
  B_tpm = c(4, 5, 18),
  D_tpm = c(7, 5, 3)
)

ggtern(emptycv,aes(A_tpm,D_tpm,B_tpm)) +
  geom_Tline(Tintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_Rline(Rintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_Lline(Lintercept = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),colour="white",linetype="dashed") +
  geom_mask() +
  labs(x="A",y="D",z="B") +
  theme_gray() +
  theme(axis.title = element_text(size = 15),axis.text = element_text(size = 10))

ggsave('fig2_panelC.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

# panel D
fig2_panelD <- read.table('fig2-panelD',header = T,sep = '\t')

heatPCsumbartab <- fig2_panelD[,1:3]
heatPCsumbartab$perc <- (heatPCsumbartab$counts/50)*100

hheatPC <- heatPCsumbartab %>%
  dplyr::select(group_id,category,perc) %>%
  pivot_wider(names_from = category,values_from = perc)
hheatPC <- as.data.frame(hheatPC)
rownames(hheatPC) <- hheatPC$group_id
hheatPC <- hheatPC[,-1]

transposed_data <- t(hheatPC)
svg('fig2_panelD.svg',width = 12,height = 9)
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


# panel E
fig2_panelE <- read.table('fig2-panelE',header = T,sep = '\t')

heatPWsumbartab <- fig2_panelE[,1:3]
heatPWsumbartab$perc <- (heatPWsumbartab$counts/50)*100

hheatPW <- heatPWsumbartab %>%
  dplyr::select(group_id,category,perc) %>%
  pivot_wider(names_from = category,values_from = perc)
hheatPW <- as.data.frame(hheatPW)
rownames(hheatPW) <- hheatPW$group_id
hheatPW <- hheatPW[,-1]

transposed_dataPW <- t(hheatPW)
svg('fig2_panelE.svg',width = 12,height = 9)
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

# panel F

fig2_panelF <- read.table('fig2-panelF',header = T,sep = '\t')
fig2_panelF$cat_over_15 <- factor(fig2_panelF$cat_over_15,levels=c('Uncategorised','DFO_b','DFO_a','DFB','Conserved'))
fig2_panelF$cross <- factor(fig2_panelF$cross,levels=c('PxC','PxW','Overlap'))

#flipped axis
ggplot(fig2_panelF) + geom_bar(aes(cat_over_15,fill=cross),colour='black',position = position_dodge2(reverse=T)) +
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
ggsave('fig2_panelF.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
