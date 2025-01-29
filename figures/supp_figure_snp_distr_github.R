library(ggplot2)

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
