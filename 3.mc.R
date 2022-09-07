#-----fig 3------
pg.mc.tpm <- tpm[,-c(8:10)]
pg.pgmc.tpm <- tpm[,c(1,11:19)]
##-----masigpro test------
BiocManager::install('maSigPro',lib='/home/nico/R/x86_64-pc-linux-gnu-library/4.0')
library(maSigPro)
data("data.abiotic")
head(data.abiotic)
data("edesign.abiotic")
head(data.abiotic)
edesign.abiotic <- as.data.frame(edesign.abiotic)

tc.test <- maSigPro (na.omit(data.abiotic), edesign.abiotic, degree = 2, vars = "groups", main = "Test")

see.genes(tc.test$sig.genes$ColdvsControl,newX11 = F,alfa = 0.05, k = 5)
STMDE66 <- data.abiotic[rownames(data.abiotic)=="STMDE66", ]
PlotGroups (STMDE66, 
            edesign = edesign.abiotic)

#----pg.mc.cor-----

## get tpm genes > 0.5 in each group 
mc.pick.genes <- unique(c(rownames(pg.pgmc.tpm[(matrixStats::rowMedians(as.matrix(pg.pgmc.tpm[,c(2:(ncol(pg.pgmc.tpm)-1))]),cols = 1:3) > 0.5) ,]), rownames(pg.pgmc.tpm[(matrixStats::rowMedians(as.matrix(pg.pgmc.tpm[,c(2:(ncol(pg.pgmc.tpm)-1))]),cols = 4:6) > 0.5) ,]), rownames(pg.pgmc.tpm[(matrixStats::rowMedians(as.matrix(pg.pgmc.tpm[,c(2:(ncol(pg.pgmc.tpm)-1))]),cols = 7:8) > 0.5) ,])))

pg.pgmc.tpm.clean <- pg.pgmc.tpm[mc.pick.genes,c(2:(ncol(pg.pgmc.tpm)-1))]
pg.pgmc.tpm.clean <- pg.pgmc.tpm.clean[,c(1:3,7,8,4:6)]

pg.mc.tpm_cor <- cor(pg.pgmc.tpm.clean,method = 'pearson',use = 'pairwise.complete.ob')
pheatmap(pg.mc.tpm_cor,border_color = NA, clustering_method = 'average',color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .9),angle_col = 315)


#-----mc.genes&plu.genes-----
msc.genes <- c('PAX7','MYOD','MYF5','MYOG','MYMK','ENSSSCG00000029441','MYH3','MYH8','DES')
diff.genes.tpm <- na.omit(tpm[c(msc.genes,plu.gene),11:18])
diff.genes.tpm <- diff.genes.tpm[,c(1:3,7,8,4:6)]
gene.anno <- data.frame(Gene_type = c(rep('Muscle differentiation',9),rep('Pluripotency',12)))
rownames(gene.anno) <- rownames(diff.genes.tpm)
anno.col <- list(Gene_type=c('Muscle differentiation'=cols[13], 'Pluripotency'=cols[5]))


pheatmap(t(scale(t(diff.genes.tpm[-c(11,nrow(diff.genes.tpm)),]))),
         border_color = NA, 
         clustering_method = 'single',
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = 315,
         annotation_row = gene.anno,
         annotation_colors = anno.col,
         breaks = unique(c(seq(-2,2, length=256))),
         main = 'TPM-filter REST_SOX2')

#-----PCA-----
group_list_pg <- factor(c( rep('pgEpiSCs.Low',3), rep('pgEpiSCs.MPC',2),rep('pgEpiSCs.MC',3)),levels = c('pgEpiSCs.Low','pgEpiSCs.MPC','pgEpiSCs.MC'))

d_p <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  # fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 70),title='Dim choose')
  fviz_pca_ind(tpm.pca,
               mean.point=F,
               # repel = T,
               label = "none",
               geom.ind = c("point",'text'),
               fill.ind = tpm.t$group_list,
               palette = rev(RColorBrewer::brewer.pal(3,"RdBu")),
               legend.title = "Groups",
               pointsize = 4,
               pointshape = 21,
               col.ind = "black",
               title = 'Samples PCA'
  ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}
d_p(pg.pgmc.tpm.clean,group_list_pg)



#----pg.mc.masigpro------

# best kmeans value 6
nk=2:15
set.seed(222)
Wss<-sapply(nk,function(k){
  kmeans(pg.pgmc.tpm.clean,centers = k,iter.max = 100)$tot.withinss})
plot(nk,Wss,type = "l",xlab="Number of k",ylab="Within sum of squares")
plot(nk,Wss,type = "o",xlab="Number of k",ylab="Within sum of squares",col='red3')


mc.design <- data.frame(Time = c(rep(0,3),rep(1,2),rep(2,3)),
                        Replicate	= c(rep(1,3),rep(2,2),rep(3,3)),
                        Group = rep(1,8))
rownames(mc.design) <- colnames(pg.pgmc.tpm.clean)


mc.masig <- maSigPro(na.omit(pg.pgmc.tpm.clean), mc.design, degree = 2, vars = "groups", main = "PG.low_pg.MC")

see.genes(mc.masig$sig.genes$Group,newX11 = F,alfa = 0.05, k = 6)

#-----pg.mc.median for STEM-----
pg.pgmc.tpm.clean.median <- data.frame(pg.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 1:3),
                                       pg.mpc.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 4:5),
                                       pg.mc.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 6:8))
rownames(pg.pgmc.tpm.clean.median) <- rownames(pg.pgmc.tpm.clean)
export(pg.pgmc.tpm.clean.median,file = '../pg.pgmc.tpm.clean.median.txt',sep = '\t',row.names = T)

#-----STEM CLUSTER--------
### cluster
## c0_p15 pass
c0_p15 <- read.table('../../MSC/c0_p15.genes.txt',header = T,sep = '\t')
c0_p15.tpm <- pg.pgmc.tpm.clean[c0_p15$symbol,]
mismatch <- c0_p15[rownames(c0_p15.tpm) %>% grep(pattern = 'NA'),] # not very important genes
mismatch.tpm <- pg.pgmc.tpm.clean[mismatch$symbol,]

c0_p15.tpm <- na.omit(pg.pgmc.tpm.clean[c0_p15$symbol,])

pheatmap(t(scale(t(c0_p15.tpm))),
         border_color = NA, 
         clustering_method = 'ward.D',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = 315,
         breaks = unique(c(seq(-2,2, length=256))),
         main = 'C0_P15.TPM')

## c0_p13 pass
c0_p13 <- read.table('../../MSC/c0_p13.genes.txt',header = T,sep = '\t')

c0_p13.tpm <- na.omit(pg.pgmc.tpm.clean[c0_p13$symbol,])

pheatmap(t(scale(t(c0_p13.tpm))),
         border_color = NA, 
         clustering_method = 'ward.D',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = 315,
         breaks = unique(c(seq(-2,2, length=256))),
         main = 'C0_P13.TPM')

## c0_p12 pass
c0_p12 <- read.table('../../MSC/c0_p12.genes.txt',header = T,sep = '\t')

c0_p12.tpm <- na.omit(pg.pgmc.tpm.clean[c0_p12$symbol,])

pheatmap(t(scale(t(c0_p12.tpm))),
         border_color = NA, 
         clustering_method = 'ward.D',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = 315,
         breaks = unique(c(seq(-2,2, length=256))),
         main = 'C0_P12.TPM')

## c0_p11 pass
c0_p11 <- read.table('../../MSC/c0_p11.genes.txt',header = T,sep = '\t')

c0_p11.tpm <- na.omit(pg.pgmc.tpm.clean[c0_p11$symbol,])

pheatmap(t(scale(t(c0_p11.tpm))),
         border_color = NA, 
         clustering_method = 'ward.D',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = 315,
         breaks = unique(c(seq(-2,2, length=256))),
         main = 'C0_P11.TPM')

## c1_p8 
c1_p8 <- read.table('../../MSC/c1_p8.genes.txt',header = T,sep = '\t')

c1_p8.tpm <- na.omit(pg.pgmc.tpm.clean[c1_p8$symbol,])

pheatmap(t(scale(t(c1_p8.tpm))),
         border_color = NA, 
         clustering_method = 'ward.D',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = 315,
         breaks = unique(c(seq(-2,2, length=256))),
         main = 'C1_P8.TPM')

## c2_p2
c2_p2 <- read.table('../../MSC/c2_p2.genes.txt',header = T,sep = '\t')

c2_p2.tpm <- na.omit(pg.pgmc.tpm.clean[c2_p2$symbol,])

pheatmap(t(scale(t(c2_p2.tpm))),
         border_color = NA, 
         clustering_method = 'ward.D',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = 315,
         # breaks = unique(c(seq(-2,2, length=256))),
         main = 'C2_P2.TPM')

## cluster integrate
cluster_inte <- rbind(c0_p15.tpm,c0_p13.tpm,c0_p12.tpm,c0_p11.tpm,c1_p8.tpm,c2_p2.tpm)

cluster_inte$Cluster <- c(rep('Cluster1',nrow(c0_p15.tpm)+nrow(c0_p13.tpm)+nrow(c0_p12.tpm)+nrow(c0_p11.tpm)),
                          rep('Cluster2',nrow(c1_p8.tpm)),
                          rep('Cluster3',nrow(c2_p2.tpm)))
cluster_inte$Profile <- c(rep('Profile_15',nrow(c0_p15.tpm)),
                          rep('Profile_13',nrow(c0_p13.tpm)),
                          rep('Profile_12',nrow(c0_p12.tpm)),
                          rep('Profile_11',nrow(c0_p11.tpm)),
                          rep('Profile_8',nrow(c1_p8.tpm)),
                          rep('Profile_2',nrow(c2_p2.tpm)))

library(ComplexHeatmap)

genes <- c('ITGB1',"MYOG","PAX7","DES","MYH3","OTX2","LIN28A","RIF1","EYA2","EVC","ZEB2","MEF2C","MYH7","RUNX1","SOX9")
gene_pos <- match(genes,rownames(cluster_inte)) %>% na.omit()
row_anno <- rowAnnotation(gene=anno_mark(at=gene_pos,labels = genes))

ComplexHeatmap::pheatmap(t(scale(t(cluster_inte[,c(1:8)]))),
         border_color = NA, 
         clustering_method = 'ward.D',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = '315',
         annotation_row = cluster_inte[,c(9:10)],
         breaks = unique(c(seq(-2,2, length=256))),
         right_annotation = row_anno)

cluster_inte_orth <- cluster_inte %>% rownames_to_column( var = 'genename') %>% left_join(y = orth[,c(1:4)], by = 'genename')
export(cluster_inte_orth,file = '../../MSC/cluster_inte_orth.tpm.xlsx') 

## STEM cluster info
library(sunburstR)

stem.info <- data.frame(v1 = c('PCG-exp.PCG-STEM.pick-sig.cluster-cluster1-profile15',
                               'PCG-exp.PCG-STEM.pick-sig.cluster-cluster1-profile11',
                               'PCG-exp.PCG-STEM.pick-sig.cluster-cluster1-profile12',
                               'PCG-exp.PCG-STEM.pick-sig.cluster-cluster1-profile13',
                               'PCG-exp.PCG-STEM.pick-sig.cluster-cluster2-profile8',
                               'PCG-exp.PCG-STEM.pick-sig.cluster-cluster3-profile2',
                               'PCG-exp.PCG-STEM.pick-unsig.cluster',
                               'PCG-exp.PCG-STEM.pick',
                               'PCG-exp.PCG',
                               'PCG-not.exp'),
                        v2 = c(nrow(c0_p15),
                               nrow(c0_p11),
                               nrow(c0_p12),
                               nrow(c0_p13),
                               nrow(c1_p8),
                               nrow(c2_p2),
                               3926,
                               8034,
                               14546,
                               nrow(data.clean)-14546
                               ))

sund2b(stem.info,
       colors = htmlwidgets::JS("d3.scaleOrdinal(d3.schemeCategory20b)")
)
sunburst(stem.info,
         count = F
)

#------


##------WGCNA for cluster P2/P11/P15/P13-P12/P8 FAILED------

library(WGCNA)
enableWGCNAThreads()
options(stringsAsFactors = FALSE)

mc.meta <- data.frame(stage = c( rep('pgEpiSCs.Low',3), rep('pgEpiSCs.MPC',2),rep('pgEpiSCs.MC',3)))
rownames(mc.meta) <- colnames(c2_p2.tpm)

## c0_p11
datExpr <- t(c0_p11.tpm) 

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# β value
powers = c(c(1:10), seq(from = 12, to=500, by=4))
#设置beta值的取值范围
sft = pickSoftThreshold(datExpr, RsquaredCut = 0.9,powerVector = powers, verbose = 5)


net = blockwiseModules(
  datExpr,
  power = 7,             #软阈值，前面计算出来的
  maxBlockSize = 6000,                   #最大block大小，将所有基因放在一个block中
  TOMType = "unsigned",                  #选择unsigned，使用标准TOM矩阵
  deepSplit = 2, minModuleSize = 30,     #剪切树参数，deepSplit取值0-4
  mergeCutHeight = 0.25,                 # 模块合并参数，越大模块越少
  numericLabels = TRUE,                  # T返回数字，F返回颜色
  pamRespectsDendro = FALSE,  
  saveTOMs = TRUE,
  saveTOMFileBase = "C0_P11-STEM-SIG-TPM-TOM",
  loadTOMs = TRUE,
  verbose = 3
)


# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
moduleColors=mergedColors

TOM = TOMsimilarityFromExpr(datExpr, power = 7);
module = "turquoise";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因
inModule = (moduleColors==module)
modProbes = probes[inModule]
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste('C0_p11', collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste('C0_p11', collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes,
  nodeAttr = moduleColors[inModule]
)

## c2

datExpr <- t(c2_p2.tpm[(matrixStats::rowMedians(as.matrix(c2_p2.tpm),cols = 1:3) > 5 & matrixStats::rowMedians(as.matrix(c2_p2.tpm),cols = 6:8) < .5) ,]) 

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# β value
powers = c(c(1:10), seq(from = 12, to=500, by=4))
#设置beta值的取值范围
sft = pickSoftThreshold(datExpr, RsquaredCut = 0.9,powerVector = powers, verbose = 5)


net = blockwiseModules(
  datExpr,
  power = 1,             #软阈值，前面计算出来的
  maxBlockSize = 6000,                   #最大block大小，将所有基因放在一个block中
  TOMType = "unsigned",                  #选择unsigned，使用标准TOM矩阵
  deepSplit = 2, minModuleSize = 30,     #剪切树参数，deepSplit取值0-4
  mergeCutHeight = 0.25,                 # 模块合并参数，越大模块越少
  numericLabels = TRUE,                  # T返回数字，F返回颜色
  pamRespectsDendro = FALSE,  
  saveTOMs = TRUE,
  saveTOMFileBase = "C0_P11-STEM-SIG-TPM-TOM",
  loadTOMs = TRUE,
  verbose = 3
)


# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
moduleColors=mergedColors

TOM = TOMsimilarityFromExpr(datExpr, power = 4);
module = "turquoise";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因
inModule = (moduleColors==module)
modProbes = probes[inModule]
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-", paste('C2_p2_8', collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-", paste('C2_p2_8', collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.5,
  nodeNames = modProbes,
  nodeAttr = moduleColors[inModule]
)


#----- -----

#-----kmeans cluster-----

# best kmeans value 6 for every samples
pg.pgmc.tpm.clean.median <- data.frame(pg.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 1:3),
                                       pg.mpc.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 4:5),
                                       pg.mc.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 6:8))
rownames(pg.pgmc.tpm.clean.median) <- rownames(pg.pgmc.tpm.clean)
pg.pgmc.tpm.clean.median.log2 <- log2(pg.pgmc.tpm.clean.median+1)

kmeansBIC = function(fit){
  
  m = ncol(fit$centers)
  n = length(fit$cluster)
  k = nrow(fit$centers)
  D = fit$tot.withinss
  return(data.frame(AIC = D + 2*m*k,
                    BIC = D + log(n)*m*k))
}

library(factoextra)
nk=2:100
set.seed(222)
Wss<-sapply(nk,function(k){
  kmeans(pg.pgmc.tpm.clean,centers = k,iter.max = 100)$tot.withinss})
# plot(nk,Wss,type = "l",xlab="Number of k",ylab="Within sum of squares")
plot(nk,Wss,type = "o",xlab="Number of k",ylab="Within sum of squares",col='red3')

# best kmeans value
Wss<-sapply(nk,function(k){
  kmeans(pg.pgmc.tpm.clean.median.log2[,c(1:3)],centers = k,iter.max = 100)$tot.withinss})
# plot(nk,Wss,type = "l",xlab="Number of k",ylab="Within sum of squares")
plot(nk,Wss,type = "o",xlab="Number of k",ylab="Within sum of squares",col='red3')


fit <- kmeans(x = pg.pgmc.tpm.clean.median.log2[,c(1:3)],10,iter.max = 100)
kmeansBIC(fit)
table(fit$cluster)
pg.pgmc.tpm.clean.median.log2$kmeans.cluster <- fit$cluster
pg.pgmc.tpm.clean.median.log2$gene <- rownames(pg.pgmc.tpm.clean.median)

kmectocluster <- data.frame(cluster = fit$cluster)

kmectocluster <- 
  kmectocluster %>% 
  rownames_to_column(var = 'gene') %>%
  left_join(y = pg.pgmc.tpm.clean.median.log2,
            by = 'gene')
kmectocluster$cluster <- paste0('C_',kmectocluster$cluster)
kmectocluster <- kmectocluster %>% group_by(cluster) %>%  mutate(cluster=paste0(cluster,"(",n()," genes)"))

kmectocluster %>% 
  select('pg.median','pg.mpc.median','pg.mc.median','cluster') %>%
  gather(key = 'stage',value = 'exp', -cluster) %>%
  mutate(stage=fct_relevel(stage, 'pg.median','pg.mpc.median','pg.mc.median')) %>%
  ggplot( aes(x=stage, y=round(exp,4),group=1,fill=cluster)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               #width=.2, 
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line") +
  labs(x="Stage", y="Expression level(log2(TPM))") +
  theme_bw() +
  theme(legend.position="none") +
  scale_fill_manual(values = rep('red',50),
                    guide=guide_legend(direction="vertical",
                                       label.position="right",
                                       title=NULL,
                                       ncol=6,
                                       label.hjust=0.8))+
  scale_color_manual(values =  rep('red',50),guide = 'none')+
  # geom_smooth()+
  # ylim(0,1)+
  facet_wrap(~cluster,scales = 'free_y',ncol = 5,drop = F)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 


### kmeans for scale data
# 
pg.pgmc.tpm.clean.median <- data.frame(pg.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 1:3),
                                       pg.mpc.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 4:5),
                                       pg.mc.median = matrixStats::rowMedians(as.matrix(pg.pgmc.tpm.clean),cols = 6:8))
rownames(pg.pgmc.tpm.clean.median) <- rownames(pg.pgmc.tpm.clean)
pg.pgmc.tpm.clean.median.scale <- data.frame(t(scale(t(pg.pgmc.tpm.clean.median))))

nk=2:100
Wss<-sapply(nk,function(k){
  kmeans(pg.pgmc.tpm.clean.median.scale,centers = k,iter.max = 100)$tot.withinss})
plot(nk,Wss,type = "o",xlab="Number of k",ylab="Within sum of squares",col='red3',)
abline(v=10,col='black')
legend("topleft", legend = 'k = 10', 
       col= "red3",
       pch = 15, bty = "n", pt.cex = 2, cex = 1.2,  horiz = F, inset =  0.1)



fit <- kmeans(x = pg.pgmc.tpm.clean.median.scale[,c(1:3)],10,iter.max = 100)
kmeansBIC(fit)
table(fit$cluster)
pg.pgmc.tpm.clean.median.scale$kmeans.cluster <- fit$cluster
pg.pgmc.tpm.clean.median.scale$gene <- rownames(pg.pgmc.tpm.clean.median)

kmectocluster <- data.frame(cluster = fit$cluster)

pg.pgmc.tpm.clean.median.log2 <- log2(pg.pgmc.tpm.clean.median+1)
pg.pgmc.tpm.clean.median.log2$gene <- rownames(pg.pgmc.tpm.clean.median.log2)

kmectocluster <- 
  kmectocluster %>% 
  rownames_to_column(var = 'gene') %>%
  left_join(y = pg.pgmc.tpm.clean.median.log2,
            by = 'gene')
kmectocluster$cluster <- paste0('C_',kmectocluster$cluster)
kmectocluster$cluster <- factor(kmectocluster$cluster,levels = c(paste0('C_',seq(1,10))))
kmectocluster <- kmectocluster %>% group_by(cluster) %>%  mutate(cluster=paste0(cluster,"(",n()," genes)"))
table(kmectocluster$cluster)
kmectocluster$cluster <- factor(kmectocluster$cluster,levels = c('C_1(1232 genes)', 'C_2(1523 genes)', 'C_3(1253 genes)', 'C_4(1333 genes)',  'C_5(1349 genes)', 'C_6(644 genes)',  'C_7(1343 genes)','C_8(1830 genes)',  'C_9(1373 genes)','C_10(2666 genes)'))

kmectocluster %>% 
  select('pg.median','pg.mpc.median','pg.mc.median','cluster') %>%
  gather(key = 'stage',value = 'exp', -cluster) %>%
  mutate(stage=fct_relevel(stage, 'pg.median','pg.mpc.median','pg.mc.median')) %>%
  ggplot( aes(x=stage, y=round(exp,4),group=1,fill=cluster)) +
  stat_summary(fun.data="mean_cl_boot", geom="ribbon",
               #width=.2, 
               alpha=I(.2)) +
  stat_summary(fun="mean", geom="line") +
  labs(x="Stage(cluster base on scale data)", y="Expression level(log2(TPM+1))") +
  theme_bw() +
  theme(legend.position="none") +
  scale_fill_manual(values = rep('red',50),
                    guide=guide_legend(direction="vertical",
                                       label.position="right",
                                       title=NULL,
                                       ncol=6,
                                       label.hjust=0.8))+
  scale_color_manual(values =  rep('red',50),guide = 'none')+
  # geom_smooth()+
  # ylim(0,1)+
  facet_wrap(~cluster,scales = 'free_y',ncol = 5,drop = F)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(labels = c('PgEpiSC','PgEpiSC.MPC','PgEpiSC.MC'))

library(rio)
export(kmectocluster,file = '../../MSC/kmeans/kmeans.cluster.xlsx')

# kmectocluster <- import('../../MSC/kmeans/kmeans.cluster.xlsx')
pg.pgmc.tpm.clean.cluster <- pg.pgmc.tpm.clean %>% rownames_to_column(var = 'gene') %>% left_join(y = kmectocluster,by = 'gene')
# pg.pgmc.tpm.clean.cluster$kmeans.cluster <- paste0('Cluster',pg.pgmc.tpm.clean.cluster$kmeans.cluster)


## 经典基因 佐证
library(ComplexHeatmap)
library(circlize)

gene_pos <- match(rownames(diff.genes.tpm),pg.pgmc.tpm.clean.cluster$gene) %>% na.omit()
row_anno <- rowAnnotation(gene=anno_mark(at=gene_pos,labels = rownames(diff.genes.tpm)))
rownames(pg.pgmc.tpm.clean.cluster) <- pg.pgmc.tpm.clean.cluster$gene

annotation_row_matrix <- data.frame(pg.pgmc.tpm.clean.cluster[,c(10)])
rownames(annotation_row_matrix) <- rownames(pg.pgmc.tpm.clean.cluster)
colnames(annotation_row_matrix) <- 'Kmeans.cluster'
annotation_col_matrix <- data.frame(Group=c(rep('pgEpiSC',3), rep('pgEpiSC.MPC',2), rep('pgEpiSC.MC',3)))
rownames(annotation_col_matrix) = colnames(pg.pgmc.tpm.clean.cluster)[2:9]

show_col(pal_npg("nrc")(10))
ha_left <- list(Kmeans.cluster=pal_npg("nrc")(10), Group = col)
names(ha_left$Kmeans.cluster) = levels(annotation_row_matrix$Kmeans.cluster)
names(ha_left$Group) = c('pgEpiSC','pgEpiSC.MPC','pgEpiSC.MC')

ComplexHeatmap::pheatmap(t(scale(t(pg.pgmc.tpm.clean.cluster[,c(2:9)]))),
         border_color = NA, 
         clustering_method = 'ward.D2',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = '315',
         annotation_row = annotation_row_matrix,
         annotation_colors = ha_left,
         annotation_col = annotation_col_matrix,
         breaks = unique(c(seq(-2,2, length=256))),
         right_annotation = row_anno)

pg.pgmc.tpm.clean.cluster <- pg.pgmc.tpm.clean.cluster %>% left_join(y = orth[,c(1:4)],by = c('gene'='genename'))
export(pg.pgmc.tpm.clean.cluster,file = '../../MSC/kmeans/pg.pgmc.tpm.clean.cluster.xlsx')

## kmeans info
kmeans.info <- data.frame(Gene.type=c('Exp.PCG','Not.exp.PCG'), count = c(14546,6734))
kmeans.info$fraction = kmeans.info$count / sum(kmeans.info$count)
kmeans.info$ymax = cumsum(kmeans.info$fraction)
kmeans.info$ymin = c(0, head(kmeans.info$ymax, n=-1))
kmeans.info$labelPosition <- (kmeans.info$ymax + kmeans.info$ymin) / 2
kmeans.info$label <- paste0(kmeans.info$Gene.type, "\n number: ", kmeans.info$count)

ggplot(kmeans.info, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Gene.type)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=5) +
  # scale_fill_manual(values = c('#e63946','#2a9d8f'),alpha = .2) +
  scale_fill_simpsons() +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")


## 根据基因选择cluster
rownames(pg.pgmc.tpm.clean.cluster) <- pg.pgmc.tpm.clean.cluster$gene
pick.genes.for.cluster <- read.table('../../MSC/kmeans/pick.genes.for.cluster.txt',header = F,sep = '\t')
pick.genes.for.cluster$gene.type <- c(rep('Pluripotency',9), rep('Early myogenesis',8), rep('Late myogenesis',20), rep('Collagen',10))
pick.genes.for.cluster.tpm <- pg.pgmc.tpm.clean.cluster %>% filter(gene %in% pick.genes.for.cluster$V1)

pick.genes.for.cluster.tpm <- separate(pick.genes.for.cluster.tpm, col = cluster, into = c('clusters','genenum'), remove = T)
pick.genes.for.cluster.tpm$cluster <- paste0(pick.genes.for.cluster.tpm$clusters,"_",pick.genes.for.cluster.tpm$genenum)
pick.genes.for.cluster.tpm <- pick.genes.for.cluster.tpm %>% left_join( y = pick.genes.for.cluster, by = c('gene' = 'V1'))
rownames(pick.genes.for.cluster.tpm) <- pick.genes.for.cluster.tpm$gene
pick.genes.for.cluster.tpm$cluster <- factor(pick.genes.for.cluster.tpm$cluster, levels = c(paste0("C_", c(1,2,3,5,7,8,10))))
pick.genes.for.cluster.tpm <- na.omit(pick.genes.for.cluster.tpm[pick.genes.for.cluster$V1,])
pick.genes.for.cluster.tpm$gene.type <- factor(pick.genes.for.cluster.tpm$gene.type, levels = c('Pluripotency','Early myogenesis','Late myogenesis','Collagen'))

ha_left <- list(cluster=pal_npg("nrc")(10)[c(1:3,5,7,8,10)], gene.type = npg.col[13:16],Group = col)
names(ha_left$cluster) = levels(pick.genes.for.cluster.tpm$cluster)
names(ha_left$gene.type) = levels(pick.genes.for.cluster.tpm$gene.type)
names(ha_left$Group) = c('pgEpiSC','pgEpiSC.MPC','pgEpiSC.MC')

ComplexHeatmap::pheatmap(t(scale(t(pick.genes.for.cluster.tpm[,c(2:9)]))),
                         border_color = NA, 
                         clustering_method = 'ward.D2',
                         show_rownames = T,
                         cluster_cols = F,cluster_rows = F,
                         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
                         angle_col = '315',
                         annotation_row = pick.genes.for.cluster.tpm[,c(18,19)],
                         annotation_colors = ha_left,
                         annotation_col = annotation_col_matrix,
                         breaks = unique(c(seq(-2,2, length=256))),fontsize_row = 8,
                         )
# Collagen genes heatmap------
collagen.genes <- import('../../MSC/kmeans/new.plot/collagen.txt',header = F)
collagen.genes.tpm <- tpm[collagen.genes$V1,]

ha_left.col <- list(Group = unique(cols)[c(1,3:6)])
names(ha_left.col$Group) <- levels(geneinfo.all$Group)[c(1,3:6)]

ComplexHeatmap::pheatmap(t(scale(t(collagen.genes.tpm[,c(11:13,17,18,14:16,2:7)]))),
                         border_color = NA, 
                         clustering_method = 'ward.D2',
                         show_rownames = T,
                         cluster_cols = F,cluster_rows = T,
                         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
                         angle_col = '315',
                         # annotation_row = pick.genes.for.cluster.tpm[,c(18,19)],
                         annotation_colors = ha_left.col,
                         annotation_col = annotation_col_matrix,
                         # breaks = unique(c(seq(-2,2, length=256))),
                         fontsize_row = 8,
)


 
#----GOBP.custome for kmeans cluster--------
mc.gobp <- import('../../MSC/kmeans/cluster.top20.GOBP.xlsx')
mc.gobp$GeneRatio <- mc.gobp$Count/mc.gobp$Background
mc.gobp.tmp <- mc.gobp[1:20,]
mc.gobp.tmp2 <- mc.gobp.tmp[order(mc.gobp.tmp$GeneRatio),]
mc.gobp.tmp2$Description <- factor(mc.gobp.tmp2$Description, levels = mc.gobp.tmp2$Description)

ggplot(data = mc.gobp.tmp2, aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = Count,color = -LogQ)) +
  scale_color_gradient(low = col[1],high = col[3], name=expression(-log[10](padj))) +
  scale_size(range  =  c(0, 6)) +
  labs( y = 'Top 20 GO BP cluster 1/7') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10,face = "bold")) 

#-----heatmap for kmeans cluster------

cluster.tpm.tmp <- pg.pgmc.tpm.clean.cluster[c(pg.pgmc.tpm.clean.cluster$cluster=='C_1(1232 genes)' | pg.pgmc.tpm.clean.cluster$cluster=='C_7(1343 genes)'),c(2:9)]
cluster.tpm.tmp <- pg.pgmc.tpm.clean.cluster[c(pg.pgmc.tpm.clean.cluster$cluster=='C_2(1523 genes)' | pg.pgmc.tpm.clean.cluster$cluster=='C_3(1253 genes)'),c(2:9)]
cluster.tpm.tmp <- pg.pgmc.tpm.clean.cluster[c(pg.pgmc.tpm.clean.cluster$cluster=='C_5(1349 genes)' | pg.pgmc.tpm.clean.cluster$cluster=='C_8(1830 genes)'),c(2:9)]
cluster.tpm.tmp <- pg.pgmc.tpm.clean.cluster[c(pg.pgmc.tpm.clean.cluster$cluster=='C_10(2666 genes)' ),c(2:9)]


pheatmap(t(scale(t(cluster.tpm.tmp))),
         border_color = NA, 
         clustering_method = 'ward.D',
         show_rownames = F,
         cluster_cols = F,cluster_rows = T,
         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
         angle_col = '315',
         # annotation_row = pick.genes.for.cluster.tpm[,c(18,19)],
         breaks = unique(c(seq(-2,2, length=256))),fontsize_row = 8,
)

# cluster 10 in vivo-------
cluster10.tpm <- tpm[pg.pgmc.tpm.clean.cluster[pg.pgmc.tpm.clean.cluster$cluster == 'C_10(2666 genes)',]$gene, ]

ComplexHeatmap::pheatmap(t(scale(t(cluster10.tpm[,c(11:13,17,18,14:16)]))),
                         border_color = NA, 
                         clustering_method = 'ward.D2',
                         show_rownames = F,
                         cluster_cols = F,cluster_rows = T,
                         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
                         angle_col = '315',
                         # annotation_row = pick.genes.for.cluster.tpm[,c(18,19)],
                         # annotation_colors = ha_left.col,
                         # annotation_col = annotation_col_matrix,
                         breaks = unique(c(seq(-2,2, length=256))),
                         fontsize_row = 8,
)

test <- cluster10.tpm[,c(11:13,14:16,5:7)]
test <- test[rowSums(test) > 0,]
ComplexHeatmap::pheatmap(t(scale(t(test))),
                         border_color = NA, 
                         clustering_method = 'ward.D2',
                         show_rownames = F,
                         cluster_cols = F,cluster_rows = T,
                         color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(5,"RdBu"),alpha=T,bias=.8)(256)),alpha = .8),
                         angle_col = '315',
                         # annotation_row = pick.genes.for.cluster.tpm[,c(18,19)],
                         # annotation_colors = ha_left.col,
                         # annotation_col = annotation_col_matrix,
                         # breaks = unique(c(seq(-2,2, length=256))),
                         fontsize_row = 8,
)

# venn of invivo mc & pg----------
## pg.low-MC
data.clean <- data.clean %>% rownames_to_column(var = 'geneid') %>% left_join( y = pig.gene.info[,c(1,2)], by = 'geneid')
rownames(data.clean) <- data.clean$genename

group = data.clean[mc.pick.genes,c(11:13,5:7)]
group_list.use <- c(rep('pgEpiSC',3),rep('MC',3))
colData = data.frame(row.names = colnames(group),group_list = group_list.use)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','MC','pgEpiSC'))

pg.low_mc <- RES %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange(padj)
pg.low_mc$change <- factor(ifelse(pg.low_mc$padj < 0.05, ifelse(pg.low_mc$log2FoldChange > 0,'MC.UP','MC.DOWN'),'NOT change'))
pg.low_mc.sig <- dplyr::filter(pg.low_mc, padj < 0.05) %>%
  dplyr::arrange(padj)

## pg.low-pg.mc
group = data.clean[mc.pick.genes,c(11:16)]
group_list.use <- c(rep('pgEpiSC',3),rep('pg.MC',3))
colData = data.frame(row.names = colnames(group),group_list = group_list.use)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pg.MC','pgEpiSC'))

pg.low_pgmc <- RES %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange(padj)
pg.low_pgmc$change <- factor(ifelse(pg.low_pgmc$padj < 0.05, ifelse(pg.low_pgmc$log2FoldChange > 0,'pgMC.UP','pgMC.DOWN'),'NOT change'))
pg.low_pgmc.sig <- dplyr::filter(pg.low_pgmc, padj < 0.05) %>%
  dplyr::arrange(padj)


## venn for LFC > 1.5 
DEGvenn = list(pgEpiSC.MC = pg.low_mc.sig[(pg.low_mc.sig$log2FoldChange > log2(1.5) | pg.low_mc.sig$log2FoldChange < -log2(1.5)),]$geneid, 
               pgEpiSC.pgMC = pg.low_pgmc.sig[(pg.low_pgmc.sig$log2FoldChange > log2(1.5) | pg.low_pgmc.sig$log2FoldChange < -log2(1.5)),]$geneid)
DEGvenn_res = Venn(DEGvenn)
plot(DEGvenn_res,doWeights = T,type="circles")

## get orth

pg.low_mc.sig <- pg.low_mc.sig %>% 
  left_join(y = pig.gene.info[,c(1:2)], by = c('geneid' = 'genename'))
pg.low_mc.sig <- pg.low_mc.sig %>% 
  left_join(y = orth[,c(1:4)], by = c('geneid.y' = 'geneid'))

pg.low_pgmc.sig <- pg.low_pgmc.sig %>% 
  left_join(y = pig.gene.info[,c(1:2)], by = c('geneid' = 'genename'))
pg.low_pgmc.sig <- pg.low_pgmc.sig %>% 
  left_join(y = orth[,c(1:4)], by = c('geneid.y' = 'geneid'))

DEG.sig <- list(pg.low_mc.sig = pg.low_mc.sig, pg.low_pgmc.sig = pg.low_pgmc.sig)
export(DEG.sig,file = '../../MSC/kmeans/pg-MC-pgMC,deg.sig.xlsx',row.names=T)


## venn for down LFC 1.5
DEGvenn_down = list(MC_down = pg.low_mc.sig[(pg.low_mc.sig$log2FoldChange < -log2(1.5) ),]$geneid, 
                    pgMC_down = pg.low_pgmc.sig[(pg.low_pgmc.sig$log2FoldChange < -log2(1.5) ),]$geneid)
DEGvenn_down_res = Venn(DEGvenn_down)
plot(DEGvenn_down_res,doWeights = F,type="circles")

## venn for up LFC 1.5
DEGvenn_up = list(MC_up = pg.low_mc.sig[(pg.low_mc.sig$log2FoldChange > log2(1.5) ),]$geneid, 
                  pgMC_up = pg.low_pgmc.sig[(pg.low_pgmc.sig$log2FoldChange > log2(1.5) ),]$geneid)
DEGvenn_up_res = Venn(DEGvenn_up)
plot(DEGvenn_up_res,doWeights = F,type="circles")

DEG.res <- list(DEGvenn_down_MC=DEGvenn_down_res@IntersectionSets[["10"]], 
                DEGvenn_down_pgMC=DEGvenn_down_res@IntersectionSets[["01"]], 
                DEGvenn_down_common=DEGvenn_down_res@IntersectionSets[["11"]], 
                DEGvenn_up_MC=DEGvenn_up_res@IntersectionSets[["10"]], 
                DEGvenn_up_pgMC=DEGvenn_up_res@IntersectionSets[["01"]], 
                DEGvenn_up_common=DEGvenn_up_res@IntersectionSets[["11"]]
)

export(DEG.res,file = '../../MSC/kmeans/new.plot/PG_MC_PGMC_DEG.res.change.xlsx',row.names=T)
export(pg.pef.tpm,file = '../pg.pef.tpm.xlsx',row.names=T)


## venn for down FC 1.5 for orth
DEGvenn_down = list(MC_down = pg.low_mc.sig[(pg.low_mc.sig$log2FoldChange < -log2(1.5) ),]$human.genename, 
                    pgMC_down = pg.low_pgmc.sig[(pg.low_pgmc.sig$log2FoldChange < -log2(1.5) ),]$human.genename)
DEGvenn_down_res = Venn(DEGvenn_down)
plot(DEGvenn_down_res,doWeights = F,type="circles")

## venn for up FC 1.5 for orth
DEGvenn_up = list(MC_up = pg.low_mc.sig[(pg.low_mc.sig$log2FoldChange > log2(1.5) ),]$human.genename, 
                  pgMC_up = pg.low_pgmc.sig[(pg.low_pgmc.sig$log2FoldChange > log2(1.5) ),]$human.genename)
DEGvenn_up_res = Venn(DEGvenn_up)
plot(DEGvenn_up_res,doWeights = F,type="circles")

DEG.res <- list(DEGvenn_down_MC=DEGvenn_down_res@IntersectionSets[["10"]], 
                DEGvenn_down_pgMC=DEGvenn_down_res@IntersectionSets[["01"]], 
                DEGvenn_down_common=DEGvenn_down_res@IntersectionSets[["11"]], 
                DEGvenn_up_MC=DEGvenn_up_res@IntersectionSets[["10"]], 
                DEGvenn_up_pgMC=DEGvenn_up_res@IntersectionSets[["01"]], 
                DEGvenn_up_common=DEGvenn_up_res@IntersectionSets[["11"]]
)

export(DEG.res,file = '../../MSC/kmeans/new.plot/PG_MC_PGMC_DEG.res.change.orth.xlsx',row.names=T)

# pheatmap(t(scale(t(pg.pef.tpm[unique(pg.low_mc.sig$geneid,pg.low_pgmc.sig$geneid),c(2:9)]))),border_color = NA, clustering_method = 'ward.D',color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),angle_col = '315',show_rownames = F)
# pheatmap(
#   t(scale(t(pg.pef.tpm[unique(c(pg.low_mc.sig[(pg.low_mc.sig$log2FoldChange < -1.5 ),]$geneid, 
#                                 pg.low_pgmc.sig[(pg.low_pgmc.sig$log2FoldChange < -1.5 ),]$geneid, 
#                                 pg.low_mc.sig[(pg.low_mc.sig$log2FoldChange > 1.5 ),]$geneid,
#                                 pg.low_pgmc.sig[(pg.low_pgmc.sig$log2FoldChange > 1.5 ),]$geneid)),c(2:9)]))),
#   border_color = NA, clustering_method = 'ward.D',color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),angle_col = '315',show_rownames = F)
# 





