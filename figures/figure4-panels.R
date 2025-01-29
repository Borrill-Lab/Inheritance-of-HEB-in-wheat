library(ggplot2)
library(ComplexUpset)
library(svglite)
library(tidyverse)

fig4_panelAupdate <- read.table('fig4-panelA_update',header = T)

ggplot(fig4_panelAupdate,aes(y = forcats::fct_rev(Cross),fill=Subgenome)) + geom_bar(position = position_dodge2(reverse = T),color='black') +
  scale_fill_manual(values = c('#579D1C','#4B1F6F','#FF950E')) +
  theme_bw() +
  labs(y='Cross',x='SNP count') +
  xlim(0,5000) +
  theme(axis.title = element_text(size = 38),axis.text = element_text(size = 32),
        legend.title = element_text(size = 32),legend.text = element_text(size = 32)) +
  annotate("text",label = paste0(dim(fig4_panelAupdate[fig4_panelAupdate$Cross=='PxC' & fig4_panelAupdate$Subgenome=='A',])[1]),x=3800,y=2.3,size = 12) +
  annotate("text",label = paste0(dim(fig4_panelAupdate[fig4_panelAupdate$Cross=='PxC' & fig4_panelAupdate$Subgenome=='B',])[1]),x=4400,y=2.0,size = 12) +
  annotate("text",label = paste0(dim(fig4_panelAupdate[fig4_panelAupdate$Cross=='PxC' & fig4_panelAupdate$Subgenome=='D',])[1]),x=1200,y=1.7,size = 12) +
  annotate("text",label = paste0(dim(fig4_panelAupdate[fig4_panelAupdate$Cross=='PxW' & fig4_panelAupdate$Subgenome=='A',])[1]),x=3500,y=1.3,size = 12) +
  annotate("text",label = paste0(dim(fig4_panelAupdate[fig4_panelAupdate$Cross=='PxW' & fig4_panelAupdate$Subgenome=='B',])[1]),x=2400,y=1.0,size = 12) +
  annotate("text",label = paste0(dim(fig4_panelAupdate[fig4_panelAupdate$Cross=='PxW' & fig4_panelAupdate$Subgenome=='D',])[1]),x=1000,y=0.7,size = 12)


ggsave(plot = last_plot(),filename = 'fig4_panelAupdate.svg',device = 'svg',width = 12,height = 6,units = 'in',dpi = 500)


fig4_panelBcis <- read.table('fig4-panelBcisupdate',header = T)
fig4_panelBtranss <- read.table('fig4-panelBtranssupdate',header = T)
fig4_panelBtransd <- read.table('fig4-panelBtransdupdate',header = T)

cispc <- upset(fig4_panelBcis,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
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
ggsave(filename = 'fig4_panelBcis.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)

transspc <- upset(fig4_panelBtranss,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
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
ggsave(filename = 'fig4_panelBtranss.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)

transspctop <- upset(fig4_panelBtranss,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
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
ggsave(filename = 'fig4_panelBtransstop.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)

transdpc <- upset(fig4_panelBtransd,c('gene_D','gene_B','gene_A','snp_D','snp_B','snp_A'),sort_sets=F, set_sizes = F,
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
ggsave(filename = 'fig4_panelBtransd.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)



fig4_panelCcis <- read.table('fig4-panelCcisupdate',header = T)
fig4_panelCtranss <- read.table('fig4-panelCtranssupdate',header = T)
fig4_panelCtransd <- read.table('fig4-panelCtransdupdate',header = T)


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
ggsave(filename = 'fig4_panelCcis.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)

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
ggsave(filename = 'fig4_panelCtranss.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)

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
ggsave(filename = 'fig4_panelCtransstop.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)

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
ggsave(filename = 'fig4_panelCtransd.svg',device = 'svg',width = 13,height = 10,units = 'in',dpi = 500)

