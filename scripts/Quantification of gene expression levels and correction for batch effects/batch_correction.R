# Marek Glombik

# This script checks the similarity of expression data in PxC and PxW and runs batch effect correction on the data.
# Note: PxW data from second batch was not further used in the analysis due to no possible correction.

library(ggplot2)
library(tidyverse)
library(edgeR)
library(sva)
library(gridExtra)

F6countsall <- read.table("/Users/glombik/work/obj1_reanalysis/F6_lines_count.tsv",header = T)
PCF6countsall <- read.table("expression_analysis_new_data/PC_F6_lines_count.tsv",header = T)
PWF6countsall <- read.table("expression_analysis_new_data/PW_F6_lines_count.tsv",header = T)

allvcfcols <- read.table('allvcfcols',header = T)
allvcfcols <- unlist(allvcfcols)
colnames(F6countsall) <- gsub("\\.", "-", colnames(F6countsall))
F6countsall <- F6countsall[,which(colnames(F6countsall) %in% allvcfcols)]

prefixold <- "old-"
colnames(F6countsall) <- paste0(prefixold, colnames(F6countsall))

prefixnew <- "new-"
colnames(PCF6countsall) <- paste0(prefixnew, colnames(PCF6countsall))
colnames(PWF6countsall) <- paste0(prefixnew, colnames(PWF6countsall))

F6countsall$genes <- rownames(F6countsall)
PCF6countsall$genes <- rownames(PCF6countsall)
PWF6countsall$genes <- rownames(PWF6countsall)

#save the names for batch effect tests
f6oldnames <- colnames(F6countsall[,1:104])
pcf6newnames <- colnames(PCF6countsall[,1:111])
pwf6newnames <- colnames(PWF6countsall[1:48])
f6newnames <- c(pcf6newnames,pwf6newnames)
f6oldbatch <- c(rep('1',times=104)) # old = 1
f6newbatch <- c(rep('2',times=159)) # new = 2
f6oldbatchdf <- data.frame(f6oldbatch,f6oldnames)
f6oldbatchdfwide <- f6oldbatchdf %>%
  pivot_wider(names_from = f6oldnames,values_from = f6oldbatch)
f6newbatchdf <- data.frame(f6newbatch,f6newnames)
f6newbatchdfwide <- f6newbatchdf %>%
  pivot_wider(names_from = f6newnames,values_from = f6newbatch)
f6allbatchdf <- cbind(f6oldbatchdfwide,f6newbatchdfwide)
colnames(f6allbatchdf) <- gsub("-", "_", colnames(f6allbatchdf))
colnames(f6allbatchdf) <- gsub("old_", "", colnames(f6allbatchdf))
colnames(f6allbatchdf) <- gsub("new_", "", colnames(f6allbatchdf))
colnames(f6allbatchdf) <- gsub("\\.", "_", colnames(f6allbatchdf))

pcallbatchdf <- f6allbatchdf %>%
  dplyr::select(starts_with('PC'))
pwallbatchdf <- f6allbatchdf %>%
  dplyr::select(starts_with('PW'))
pcnat_order <- order(colnames(pcallbatchdf))
pcallbatchdf <- pcallbatchdf[,pcnat_order]
pwnat_order <- order(colnames(pwallbatchdf))
pwallbatchdf <- pwallbatchdf[,pwnat_order]
pcallbatchdf$id <- 'batch'
pcallbatchdf <- pcallbatchdf %>%
  dplyr::select(id,PC_1_3:PC_P_3)
pwallbatchdf$id <- 'batch'
pwallbatchdf <- pwallbatchdf %>%
  dplyr::select(id,PW_10_2:PW_W_3)

# write.table(f6allbatchdf,'f6allbatchdf',row.names = F,quote = F,sep = '\t')
# write.table(pwallbatchdf,'pwallbatchdf',row.names = F,quote = F,sep = '\t')
# write.table(pcallbatchdf,'pcallbatchdf',row.names = F,quote = F,sep = '\t')


# Continue
allF6HC <- F6countsall %>%
  inner_join(PCF6countsall) %>%
  inner_join(PWF6countsall) %>%
  dplyr::select(!genes)
rownames(allF6HC) <- F6countsall$genes

allF6HC <- allF6HC[!grepl("LC$", rownames(allF6HC)), ]
colnames(allF6HC) <- gsub("\\.", "_", colnames(allF6HC))
colnames(allF6HC) <- gsub("\\-", "_", colnames(allF6HC))
colnames(allF6HC) <- gsub("old_", "old-", colnames(allF6HC))
colnames(allF6HC) <- gsub("new_", "new-", colnames(allF6HC))

y <- allF6HC
group <- factor(c(rep(1:(length(colnames(y))),each=1)))
y <- DGEList(counts=y,group = group)
keep <- filterByExpr(y,min.count =10,group = group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method = 'TMM')
cpmF6countsallHC <- as.data.frame(cpm(y,log = F,normalized.lib.sizes = T))
cpmF6countsallHC$rowmeans <- rowMeans(cpmF6countsallHC)
cpmF6countsallHC <- cpmF6countsallHC[cpmF6countsallHC$rowmeans>1,1:263]



zcpm <- scale(cpmF6countsallHC,center = T,scale = T)
tcpm <- t(zcpm)
cpm_pca <- prcomp(tcpm)
cpm_out <- as.data.frame(cpm_pca$x)
cpm_out$group <- sapply(strsplit(as.character(row.names(tcpm)),"_"),"[[",1)
head(cpm_out)
percentage <- round(cpm_pca$sdev / sum(cpm_pca$sdev) * 100, 2)
percentage <- paste( colnames(cpm_out), "(", paste( as.character(percentage), "%", ")", sep="") )


p <- ggplot(cpm_out,aes(x=PC1,y=PC2,color=group)) +
  geom_point(size=4.5,alpha=0.5) +
  # geom_text(aes(label = rownames(cpm_out))) +
  xlab(percentage[1]) +
  ylab(percentage[2])
p

# PC only
PCcpmF6countsallHC <- cpmF6countsallHC %>%
  dplyr::select(matches('PC'))
zcpm <- scale(PCcpmF6countsallHC,center = T,scale = T)
tcpm <- t(zcpm)
cpm_pca <- prcomp(tcpm)
cpm_out <- as.data.frame(cpm_pca$x)
cpm_out$group <- sapply(strsplit(as.character(row.names(tcpm)),"_"),"[[",1)
head(cpm_out)
percentage <- round(cpm_pca$sdev / sum(cpm_pca$sdev) * 100, 2)
percentage <- paste( colnames(cpm_out), "(", paste( as.character(percentage), "%", ")", sep="") )


q <- ggplot(cpm_out,aes(x=PC1,y=PC2,color=group)) +
  geom_point(size=4.5,alpha=0.5) +
  # geom_text(aes(label = rownames(cpm_out))) +
  xlab(percentage[1]) +
  ylab(percentage[2])



# PW only
PWcpmF6countsallHC <- cpmF6countsallHC %>%
  dplyr::select(matches('PW'))
zcpm <- scale(PWcpmF6countsallHC,center = T,scale = T)
tcpm <- t(zcpm)
cpm_pca <- prcomp(tcpm)
cpm_out <- as.data.frame(cpm_pca$x)
cpm_out$group <- sapply(strsplit(as.character(row.names(tcpm)),"_"),"[[",1)
head(cpm_out)
percentage <- round(cpm_pca$sdev / sum(cpm_pca$sdev) * 100, 2)
percentage <- paste( colnames(cpm_out), "(", paste( as.character(percentage), "%", ")", sep="") )


r <- ggplot(cpm_out,aes(x=PC1,y=PC2,color=group)) +
  geom_point(size=4.5,alpha=0.5) +
  # geom_text(aes(label = rownames(cpm_out))) +
  xlab(percentage[1]) +
  ylab(percentage[2])


grid.arrange(p,q,r,nrow=2)


write.table(cpmF6countsallHC,file = 'cpmF6countsallHC',row.names = T,quote = F,sep = '\t')
write.table(allF6HC,file = 'allF6HC',row.names = T,quote = F,sep = '\t')


allF6HC <- read.table('allF6HC',header = T)
PCallF6HC <- allF6HC %>%
  dplyr::select(matches('PC'))
PWallF6HC <- allF6HC %>%
  dplyr::select(matches('PW'))
colnames(PCallF6HC) <- gsub("old\\-", "", colnames(PCallF6HC))
colnames(PCallF6HC) <- gsub("new\\-", "", colnames(PCallF6HC))
colnames(PWallF6HC) <- gsub("old\\-", "", colnames(PWallF6HC))
colnames(PWallF6HC) <- gsub("new\\-", "", colnames(PWallF6HC))
# write.table(PCallF6HC,file = 'PCallF6HC',row.names = T,quote = F,sep = '\t')
# write.table(PWallF6HC,file = 'PWallF6HC',row.names = T,quote = F,sep = '\t')

cpmF6countsallHC <- read.table('cpmF6countsallHC',header = T)
PCallF6HC <- cpmF6countsallHC %>%
  dplyr::select(matches('PC'))
PWallF6HC <- cpmF6countsallHC %>%
  dplyr::select(matches('PW'))
colnames(PCallF6HC) <- gsub("old\\.", "", colnames(PCallF6HC))
colnames(PCallF6HC) <- gsub("new\\.", "", colnames(PCallF6HC))
colnames(PWallF6HC) <- gsub("old\\.", "", colnames(PWallF6HC))
colnames(PWallF6HC) <- gsub("new\\.", "", colnames(PWallF6HC))
# write.table(PCallF6HC,file = 'PCcpmF6HC',row.names = T,quote = F,sep = '\t')
# write.table(PWallF6HC,file = 'PWcpmF6HC',row.names = T,quote = F,sep = '\t')

# Correct for batch effect with Combat-seq

# load in the batch data info

f6allbatchdf <- read.table('f6allbatchdf',header=T)
f6allbatchdf <- as.data.frame(t(f6allbatchdf))
colnames(f6allbatchdf) <- c('batch')
f6allbatchdf$cross <- rownames(f6allbatchdf)
f6allbatchdf$cross <- substr(f6allbatchdf$cross,1,2)
f6allbatchdf$cross[f6allbatchdf$cross == "PC"] <- 1
f6allbatchdf$cross[f6allbatchdf$cross == "PW"] <- 2

#adjust the expression data colnames
colnames(allF6HC) <- gsub("old.", "", colnames(allF6HC))
colnames(allF6HC) <- gsub("new.", "", colnames(allF6HC))

f6batch = f6allbatchdf$batch
f6cross = as.integer(f6allbatchdf$cross)

### Now run Combat-seq 
f6combseqmat <- as.matrix(allF6HC)
f6combatseq_edata = ComBat_seq(f6combseqmat,batch = f6batch,group = f6cross)
f6combatseq_edata <- as.data.frame(f6combatseq_edata)
head(f6combatseq_edata)

y <- f6combatseq_edata
group <- factor(c(rep(1:(length(colnames(y))),each=1)))
y <- DGEList(counts=y,group = group)
keep <- filterByExpr(y,min.count =10,group = group)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y,method = 'TMM')
cpmF6countsallHC <- as.data.frame(cpm(y,log = F,normalized.lib.sizes = T))
cpmF6countsallHC$rowmeans <- rowMeans(cpmF6countsallHC)
cpmF6countsallHC <- cpmF6countsallHC[cpmF6countsallHC$rowmeans>1,1:263]
write.table(cpmF6countsallHC,file='cpmf6combatseq_edata',quote=F,sep='\t',row.names=T)
