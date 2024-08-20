# This script produces panels for supplementary figure 1 from intermediate data tables created during the analysis
library(ggplot2)
library(svglite)

# Panel A
suppfig1_panelA <- read.table('/Users/glombik/work/Documents/my_articles/inheritance_2023/figure_scripts_github/suppfig1-panelA',header = T)

ggplot(suppfig1_panelA) +
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

ggsave(filename = 'suppfig1_panelA.svg',device = 'svg',dpi = 400,width = 14,height = 10,units = "in")

# Panel B
suppfig1_panelB <- read.table('/Users/glombik/work/Documents/my_articles/inheritance_2023/figure_scripts_github/suppfig1-panelB',header = T)

ggplot(suppfig1_panelB,aes(cross.x,distance,color=cross.x)) +
  geom_violin(show.legend = F) +
  geom_point(position = 'jitter',size=1,alpha=0.6,show.legend = F) +
  scale_color_manual(values = c('black','#009128ff')) +
  theme_classic() +
  xlab('Cross') +
  ylab('Bias distance between\ntriads in parental lines') +
  theme(axis.text = element_text(size = 32),
        axis.title = element_text(size=38))
ggsave('suppfig1_panelB.svg',device = 'svg',height = 8,width = 12,dpi = 400,units = 'in',plot = last_plot())
