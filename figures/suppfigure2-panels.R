# This script produces panels for supplementary figure 1 from intermediate data tables created during the analysis
library(ggplot2)
library(ggtern)

# Panel A
suppfig2_panelA_PC <- read.table('suppfig2-panelA-PC',header = T)
suppfig2_panelA_PW <- read.table('suppfig2-panelA-PW',header = T)

ggtern(suppfig2_panelA_PC,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='black') +
  scale_fill_gradient2(low = 'white',high = 'black') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('suppfig2_panelA_PC.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

ggtern(suppfig2_panelA_PW,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='#009128ff') +
  scale_fill_gradient2(low = 'white',high = '#009128ff') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('suppfig2_panelA_PW.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

# Panel B
suppfig2_panelB_PC <- read.table('suppfig2-panelB-PC',header = T)
suppfig2_panelB_PW <- read.table('suppfig2-panelB-PW',header = T)

ggtern(suppfig2_panelB_PC,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='black') +
  scale_fill_gradient2(low = 'white',high = 'black') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('suppfig2_panelB_PC.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

ggtern(suppfig2_panelB_PW,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='#009128ff') +
  scale_fill_gradient2(low = 'white',high = '#009128ff') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('suppfig2_panelB_PW.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())


# Panel C
suppfig2_panelC_PC <- read.table('suppfig2-panelC-PC',header = T)
suppfig2_panelC_PW <- read.table('suppfig2-panelC-PW',header = T)

ggtern(suppfig2_panelC_PC,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='black') +
  scale_fill_gradient2(low = 'white',high = 'black') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('suppfig2_panelC_PC.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

ggtern(suppfig2_panelC_PW,aes(A,D,B)) + 
  stat_density_tern(geom = 'polygon',aes(fill=after_stat(log(level))),bins=20,bdl = 0.01,bdl.val = NA) +
  geom_point(size=0.1,alpha=0.2,color='#009128ff') +
  scale_fill_gradient2(low = 'white',high = '#009128ff') +
  theme_void() +
  theme_nolabels() +
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.title = element_blank())
ggsave('suppfig2_panelC_PW.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())

