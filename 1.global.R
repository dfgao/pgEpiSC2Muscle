
.libPaths('/home/nico/R/x86_64-pc-linux-gnu-library/4.0')

#----load packages-------

library(DESeq2)
library(ggplot2)
library(factoextra)
library(FactoMineR)
library(reshape2)
library(ggplot2)
library(Vennerable)
library(rio)
library(ggsci)
library(RColorBrewer)
library(viridis)
library(tidyverse)
library(hrbrthemes)
library(pheatmap)
col <- c('#2a9d8f','#e9c46a','#e63946')
npg.col <- c("#E7959B","#DB5E67","#CFD99F",'#C5EC41','#E5CD94',"#FC0C00","#7C4374","#339E55","#000376","#2A82C6","#8C6D36","#CB70C0",
             "#EBB854",'#FC8D37',"#63753A","#6D2811","#DD9AD2","#68AADE","#3B397C","#9D9AE5","#B8CF6E","#949494","#BF4F8F","#844346")
options(scipen = 6)

#----load data----

### load gene info
# getwd() "/home/nico/project/01.zhugaoxiang/04.pig.RNA.MSC/05.analysis/pgEpiSC/zgx_pg"
pig.gene.info <- import('~/project/99.gaodengfeng/01.multi.embory/06.analysis/00.datainfo/01.pig/01.mydata/sus.105.xlsx',which = 'pig_105_clean')
pig.gene.info <- pig.gene.info[pig.gene.info$biotype == 'protein_coding' | pig.gene.info$biotype == 'lncRNA' | pig.gene.info$biotype == 'pseudogene' | pig.gene.info$biotype == 'processed_pseudogene',]
which(pig.gene.info$genename =='DEFB1')
pig.gene.info[10314,2] <- 'AMELY.X'
pig.gene.info[14647,2] <- 'AMELY.Y'
pig.gene.info[12197,2] <- 'GZMA.1'
pig.gene.info[12198,2] <- 'GZMA.2'
pig.gene.info = pig.gene.info[-c(14831),]

### input raw count
samples <- read.table('../../../01.data/17rnaseq/rna.sample.txt',sep = '\t')
data <- read.table('../../../04.quant/MSC.ref105.txt',header = T,sep = '\t',comment.char = '#',row.names = 1)
data <- data[pig.gene.info[pig.gene.info$biotype=='protein_coding',]$geneid, ]
data <- na.omit(data)
data.clean <- data[,c(6:ncol(data))]
colnames(data.clean) <- samples$V1
data.clean <- na.omit(data.clean)
meta <- data[,c(1:5)]

### calculate TPM 
kb <- meta$Length / 1000
rpk <- data.clean / kb
tpm <- data.frame(t(t(rpk)/colSums(rpk) * 1000000))
rm(kb,rpk)

### align ratio plot
align <- import('../../rna align ratio.xlsx')
colnames(align) <- c('Sample','Reads','Ratio')
align$Group <- factor(c( rep('pgEpiSCs.Low',3),rep('pgEpiSCs.High',3),rep('pgEpiSCs.MC',3),rep('pgEpiSCs.MPC',2),rep('MPC',3),rep('MPC.Dif',3)),levels = c('pgEpiSCs.Low','pgEpiSCs.High','pgEpiSCs.MPC','pgEpiSCs.MC','MPC','MPC.Dif'))
align <- align[c(1:6,10,11,7:9,12:17),]
align$Sample <- factor(align$Sample, levels = align$Sample)
align$reads.scale <- (round(align$Reads/10^6,1))
align$pct.scale <- align$Ratio * 25

ggplot(data = align) +
  theme_bw() +
  # theme_ipsum()+
 geom_bar(stat = "identity",aes(x=Sample,y=reads.scale,fill=Group),alpha=.8) +
  geom_line(aes(x = Sample, y = pct.scale, group = 1), size = 1, color = '#009688') +
  geom_point(aes(x = Sample, y = pct.scale, group = 1), size = 2, shape = 19, color='#762a83') +
  scale_fill_brewer(palette = "RdBu",direction = -1,) +
  labs(x=NULL) +
  scale_y_continuous(name = "Reads number (MB)",
                     limits = c(0,30),
                     breaks = seq(0, 30, 5),
                     expand = c(0, 0),
                     sec.axis = sec_axis(~./1 ,
                                         name = "Align percent (%)",
                                         breaks = seq(0, 30, 5),
                                         labels = seq(0,120,20))
  ) +
  theme(axis.text.x = element_text(angle = 45,margin = 
                                     margin(1,0,0,0,'cm')),
        axis.text.y = element_text(size = 12)
        )


#----correlation-----

library(corrplot)
tpm_cor <- cor(tpm,method = 'pearson',use = 'pairwise.complete.ob')
pheatmap(tpm_cor,border_color = NA, clustering_method = 'average',color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),angle_col = 315)

#----global PCA-------

group_list <- factor(c( rep('MPC',3),rep('MPC.Dif',3), rep('pgEpiSCs.High',3), rep('pgEpiSCs.Low',3), rep('pgEpiSCs.MC',3),rep('pgEpiSCs.MPC',2)),levels = c('pgEpiSCs.Low','pgEpiSCs.High','pgEpiSCs.MPC','pgEpiSCs.MC','MPC','MPC.Dif'))

d_p <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 40),title='Dim choose')
  fviz_pca_ind(tpm.pca,
               mean.point=F,
               # repel = T,
               label = "none",
               geom.ind = c("point",'text'),
               fill.ind = tpm.t$group_list,
               palette = rev(RColorBrewer::brewer.pal(6,"RdBu")),
               legend.title = "Groups",
               pointsize = 4,
               pointshape = 21,
               col.ind = "black",
               title = 'Samples PCA'
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}
d_p(tpm,group_list)

rm(tpm.t,var,tpm.pca)


