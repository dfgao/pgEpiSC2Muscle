#----Base on TPM filter-----
options(scipen = 10)
clean_tpm <- tpm %>% select(one_of(samples$V1)) %>% drop_na()

clean.genes <- unique(c(rownames(clean_tpm[(matrixStats::rowMedians(as.matrix(clean_tpm),cols = 1:3) > 0.5) ,]), 
                        rownames(clean_tpm[(matrixStats::rowMedians(as.matrix(clean_tpm),cols = 4:6) > 0.5) ,]), 
                        rownames(clean_tpm[(matrixStats::rowMedians(as.matrix(clean_tpm),cols = 7:9) > 0.5) ,]),
                        rownames(clean_tpm[(matrixStats::rowMedians(as.matrix(clean_tpm),cols = 10:12) > 0.5) ,]),
                        rownames(clean_tpm[(matrixStats::rowMedians(as.matrix(clean_tpm),cols = 13:15) > 0.5) ,]),
                        rownames(clean_tpm[(matrixStats::rowMedians(as.matrix(clean_tpm),cols = 16:17) > 0.5) ,])
                        ))
clean_tpm <- clean_tpm[clean.genes,] %>% drop_na() 
clean_tpm <- round(clean_tpm,3)


#----correlation-----

library(corrplot)
tpm_cor <- cor(clean_tpm,method = 'pearson',use = 'pairwise.complete.ob')
pheatmap(tpm_cor,border_color = NA, clustering_method = 'average',color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(7,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),angle_col = 315)

#----global PCA-------

group_list <- factor(c( rep('MPC',3),rep('MPC.Dif',3), rep('pgEpiSCs.High',3), rep('pgEpiSCs.Low',3), rep('pgEpiSCs.MC',3),rep('pgEpiSCs.MPC',2)),levels = c('pgEpiSCs.Low','pgEpiSCs.High','pgEpiSCs.MPC','pgEpiSCs.MC','MPC','MPC.Dif'))

d_p <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  # fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 45),title='Dim choose')
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
d_p(clean_tpm,group_list)

#-----hclsut-----
hc_p <- function(data){
  dist.r = dist(t(data), method="euclidean") 
  hc.r = hclust(dist.r, method = "ward.D") 
  plot(hc.r,main='hclust')
}
hc_p(clean_tpm)

#-----gene info-----
geneinfo <- data.frame(TPM_100 = apply(clean_tpm,2,function(x) {table(x>=100)["TRUE"]}),
                       TPM_100_50 = apply(clean_tpm,2,function(x) {table(x>=50 & x< 100)["TRUE"]}),
                       TPM_50_5 = apply(clean_tpm,2,function(x) {table(x>=5 & x< 50)["TRUE"]}),
                       TPM_5_0.5 = apply(clean_tpm,2,function(x) {table(x>=0.5 & x< 5)["TRUE"]}),
                       TPM_0.5 = apply(clean_tpm,2,function(x) {table(x< 0.5)["TRUE"]})
)
geneinfo.all <- geneinfo %>% rownames_to_column( var = 'Sample') %>% left_join(y = align,by = 'Sample')

test <- melt(geneinfo.all[,c(1:6,9)],id.vars = c('Sample','Group'),variable.name = 'Range',value.name = 'genenum')
test$Sample <- factor(test$Sample, levels = c(paste0('pgEpiSCs.Low.',seq(1,3)),
                                              paste0('pgEpiSCs.High.',seq(1,3)),
                                              paste0('pgEpiSCs.MPC.',seq(1,2)),
                                              paste0('pgEpiSCs.MC.',seq(1,3)),
                                              paste0('MPC.',seq(1,3)),
                                              paste0('MPC.Dif.',seq(1,3))
                                              ))


ggplot(data = test) +
  theme_bw() +
  # theme_ipsum()+
  geom_bar(stat = "identity",aes(x=Sample,y=genenum,fill=Range),alpha=.8,position = 'stack') +
  geom_vline(xintercept = c(3.5,6.5,8.5,11.5,14.5)) +
  scale_fill_brewer(palette = "RdBu",direction = 1,) +
  labs(y='Gene number',x=NULL) +
  theme(axis.text.x = element_text(angle = 45,margin = 
                                     margin(1,0,0,0,'cm'),size = 8),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 15)
  )

#-----gene complex-----
library(GGally)
library(grid)
clean.tpm.sort <- data.frame(apply(clean_tpm, 2, function(x){sort(x,decreasing = T)}))
# gene.comp <- data.frame(top20 = colSums(clean.tpm.sort[c(1:3006),])/colSums(clean.tpm.sort),
#                         top40 = colSums(clean.tpm.sort[c(3007:6012),])/colSums(clean.tpm.sort),
#                         top60 = colSums(clean.tpm.sort[c(round(ncol(clean.tpm.sort)*0.4,0):round(ncol(clean.tpm.sort)*0.,0)),])/colSums(clean.tpm.sort),)
gene.comp <- clean.tpm.sort

for (i in samples$V1) {
  gene.comp[,i] <- clean.tpm.sort[,i]/colSums(clean.tpm.sort)[i]
}
gene.comp <- data.frame(apply(gene.comp, 2, cumsum))

plot.data <- data.frame(t(gene.comp),Group=group_list) %>% rownames_to_column(var = 'Samples')
colnames(plot.data) <- c('Samples',(seq(1:15031)),'Group')
plot.data.log <- melt(plot.data,id.vars = c('Group','Samples'),variable.name = 'Genenum',value.name = 'Percent')
plot.data.log$Genenum <- as.numeric(plot.data.log$Genenum)
plot.data.log$Samples <- factor(plot.data.log$Sample, levels = c(paste0('pgEpiSCs.Low.',seq(1,3)),
                                                        paste0('pgEpiSCs.High.',seq(1,3)),
                                                        paste0('pgEpiSCs.MPC.',seq(1,2)),
                                                        paste0('pgEpiSCs.MC.',seq(1,3)),
                                                        paste0('MPC.',seq(1,3)),
                                                        paste0('MPC.Dif.',seq(1,3))
))
plot.data.log$Percent <- plot.data.log$Percent*100

cols <- rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=1)(6))
cols <- as.matrix(apply(data.frame(cols),1,function(x){rep(x,3)}))
cols <- as.vector(cols)
ggplot(data = plot.data.log) +
  theme_bw() +
  geom_line(aes(x=Genenum, y=Percent,color=Samples),size=.7 ) + 
  scale_x_log10() +
  scale_color_manual(values = cols[c(1:11,13:18)]) +
  labs(x = 'Genes number',y = 'Gene accumulation ratio (%)') +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.key.height=unit(.85,"line"))



