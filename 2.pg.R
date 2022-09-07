### Fig 1
library(ComplexHeatmap)
orth <- import('../../pigref105_human.unique.xlsx',which = 'Sheet1')

#----load data-------

## PEF 
pef <- read.table('../../../04.quant/PEF.ref105.txt',header = T,sep = '\t',comment.char = '#',row.names = 1)
pef <- pef[pig.gene.info[pig.gene.info$biotype=='protein_coding',]$geneid, ]
pef <- na.omit(pef)
pef.clean <- pef[,c(6:ncol(pef))]
colnames(pef.clean) <- c('PEF.1','PEF.2')
pef.clean <- na.omit(pef.clean)

## calculate TPM 
kb <- meta$Length / 1000
rpk <- pef.clean / kb
pef.tpm <- data.frame(t(t(rpk)/colSums(rpk) * 1000000))
rm(kb,rpk)

## get data for analysis
pg.pef.count <- cbind(data.clean[,c(10:12,7:9)],pef.clean)
pg.pef.tpm <- cbind(tpm[,c(10:12,7:9)],pef.tpm)


#-----A-Pluripotent genes------
pg.pef.count <- pg.pef.count %>% rownames_to_column(var = 'geneid') %>% left_join( y = pig.gene.info[,c(1,2)], by = 'geneid')
rownames(pg.pef.count) <- pg.pef.count$genename
pg.pef.tpm <- pg.pef.tpm %>% rownames_to_column(var = 'geneid') %>% left_join( y = pig.gene.info[,c(1,2)], by = 'geneid')
rownames(pg.pef.tpm) <- pg.pef.tpm$genename

plu.gene <- c('NANOG','SOX2','ENSSSCG00000001393','OTX2','LIN28A','TCF3','FGF2','LEFTY2','SMARCAD1','MYST3','SETDB1','JARID2','RIF1','REST')
tpm.plu <- na.omit(pg.pef.tpm[plu.gene,-c(1,ncol(pg.pef.tpm))])
pheatmap(t(scale(t(tpm.plu))),border_color = NA, clustering_method = 'ward.D2',color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=.8)(256)),alpha = 1),angle_col = 315,breaks = unique(c(seq(-2,2, length=256))))
pheatmap(t(scale(t(tpm.plu))),border_color = NA, cluster_rows = T,cluster_cols = F,color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=.8)(256)),alpha = 1),angle_col = '315',breaks = unique(c(seq(-2,2, length=256))))



# tpm.plu.bar <- log2(tpm.plu+1)
# tpm.plu.bar %>% 
#   gather(key = 'Samples',value="Val") %>% 
#   mutate(Groups = factor(c(rep('pgEpiSCs.Low',(nrow(tpm.plu)+1)*3 ), rep('pgEpiSCs.High',(nrow(tpm.plu)+1)*3 ), rep('PEF',(nrow(tpm.plu)+1)*2 )),levels = c('pgEpiSCs.Low','pgEpiSCs.High','PEF'))) %>% 
#   ggplot( aes(fill=Groups, y=Val, x=Samples)) + 
#   geom_violin() +
#   geom_jitter(color="black", size=0.4, alpha=0.9)

# tpm.plu.bar <- log10(tpm.plu+1)
# tpm.plu.bar %>% 
#   gather(key = 'Samples',value="Val") %>% 
#   mutate(Groups = factor(c(rep('pgEpiSCs.Low',(nrow(tpm.plu))*3 ), rep('pgEpiSCs.High',(nrow(tpm.plu))*3 ), rep('PEF',(nrow(tpm.plu))*2 )),levels = c('pgEpiSCs.Low','pgEpiSCs.High','PEF'))) %>% 
#   ggplot( aes(fill=Groups, y=Val, x=Samples)) + 
#   geom_violin() +
#   geom_jitter(color="black", size=0.4, alpha=0.9)+
#   theme_ipsum()


#-----S1B_correlation------
pg.tpm_cor <- cor(pg.pef.tpm[,c(2:9)],method = 'pearson',use = 'pairwise.complete.ob')
pheatmap(pg.tpm_cor,border_color = NA, clustering_method = 'average',color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(8,"RdBu"),alpha=T,bias=1)(256)),alpha = .9),angle_col = 315)


## get tpm genes > 0.5 in each group 
pick.genes <- unique(c(rownames(pg.pef.tpm[(matrixStats::rowMedians(as.matrix(pg.pef.tpm[,c(2:9)]),cols = 1:3) > 0.5) ,]), rownames(pg.pef.tpm[(matrixStats::rowMedians(as.matrix(pg.pef.tpm[,c(2:9)]),cols = 4:6) > 0.5) ,]), rownames(pg.pef.tpm[(matrixStats::rowMedians(as.matrix(pg.pef.tpm[,c(2:9)]),cols = 7:8) > 0.5) ,])))


#-----S1A_PCA------
group_list_pg <- factor(c( rep('pgEpiSCs.Low',3), rep('pgEpiSCs.High',3),rep('PEF',2)),levels = c('pgEpiSCs.Low','pgEpiSCs.High','PEF'))

d_p <- function(tpm,group_list){
  tpm.t = as.data.frame(t(tpm))
  tpm.t = cbind(tpm.t,group_list)
  tpm.pca <- PCA(tpm.t[,-ncol(tpm.t)],graph = FALSE)
  fviz_screeplot(tpm.pca, addlabels = TRUE, ylim = c(0, 70),title='Dim choose')
  # fviz_pca_ind(tpm.pca,
  #              mean.point=F,
  #              # repel = T,
  #              label = "none",
  #              geom.ind = c("point",'text'),
  #              fill.ind = tpm.t$group_list,
  #              palette = rev(RColorBrewer::brewer.pal(3,"RdBu")),
  #              legend.title = "Groups",
  #              pointsize = 4,
  #              pointshape = 21,
  #              col.ind = "black",
  #              title = 'Samples PCA'
  # ) + theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15))
}
d_p(pg.pef.tpm[pick.genes,c(2:9)],group_list_pg)


#-----B-DEG------- 
## pg.low-PEF
group = pg.pef.count[pick.genes,c(2:4,8:9)]
group_list.use <- c(rep('pgEpiSCs.Low',3),rep('PEF',2))
colData = data.frame(row.names = colnames(group),group_list = group_list.use)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pgEpiSCs.Low','PEF'))

pg.low_pef <- RES %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange(padj)
pg.low_pef$change <- factor(ifelse(pg.low_pef$padj < 0.05, ifelse(pg.low_pef$log2FoldChange > 0,'pg.Low.UP','pg.Low.DOWN'),'NOT change'))
pg.low_pef.sig <- dplyr::filter(pg.low_pef, padj < 0.05) %>%
  dplyr::arrange(padj)

## pg.high-PEF
group = pg.pef.count[pick.genes,c(5:9)]
group_list.use <- c(rep('pgEpiSCs.High',3),rep('PEF',2))
colData = data.frame(row.names = colnames(group),group_list = group_list.use)

dds <- DESeqDataSetFromMatrix(countData = group,colData = colData,design = ~group_list)
dds <- DESeq(dds)
plotDispEsts(dds)
RES <- results(dds, contrast = c('group_list','pgEpiSCs.High','PEF'))

pg.high_pef <- RES %>%
  data.frame() %>%
  rownames_to_column(var="geneid") %>%
  as_tibble() %>%
  arrange(padj)
pg.high_pef$change <- factor(ifelse(pg.high_pef$padj < 0.05, ifelse(pg.high_pef$log2FoldChange > 0,'pg.High.UP','pg.High.DOWN'),'NOT change'))
pg.high_pef.sig <- dplyr::filter(pg.high_pef, padj < 0.05) %>%
  dplyr::arrange(padj)


## venn for LFC > 1.5 
DEGvenn = list(pgEpiSC.Low_pEF = pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange > 1.5 | pg.low_pef.sig$log2FoldChange < -1.5),]$geneid, 
               pgEpiSC.High_pEF = pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange > 1.5 | pg.high_pef.sig$log2FoldChange < -1.5),]$geneid)
DEGvenn_res = Venn(DEGvenn)
plot(DEGvenn_res,doWeights = F,type="circles")

## get orth

pg.low_pef.sig <- pg.low_pef.sig %>% 
  left_join(y = pig.gene.info[,c(1:2)], by = c('geneid' = 'genename'))
pg.low_pef.sig <- pg.low_pef.sig %>% 
  left_join(y = orth[,c(1:4)], by = c('geneid.y' = 'geneid'))

pg.high_pef.sig <- pg.high_pef.sig %>% 
  left_join(y = pig.gene.info[,c(1:2)], by = c('geneid' = 'genename'))
pg.high_pef.sig <- pg.high_pef.sig %>% 
  left_join(y = orth[,c(1:4)], by = c('geneid.y' = 'geneid'))

DEG.sig <- list(pg.low_pef.sig = pg.low_pef.sig, pg.high_pef.sig = pg.high_pef.sig)
export(DEG.sig,file = '../DEG.sig.xlsx',row.names=T)



## venn for down FC 1.5
DEGvenn_down = list(pgEpiSC.Low_down = pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange < -1.5 ),]$geneid, 
               pgEpiSC.High_down = pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange < -1.5 ),]$geneid)
DEGvenn_down_res = Venn(DEGvenn_down)
plot(DEGvenn_down_res,doWeights = F,type="circles")

## venn for up FC 1.5
DEGvenn_up = list(pgEpiSC.Low_up = pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange > 1.5 ),]$geneid, 
                    pgEpiSC.High_up = pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange > 1.5 ),]$geneid)
DEGvenn_up_res = Venn(DEGvenn_up)
plot(DEGvenn_up_res,doWeights = F,type="circles")

DEG.res <- list(DEGvenn_down_pg.low=DEGvenn_down_res@IntersectionSets[["10"]], 
                DEGvenn_down_pg.high=DEGvenn_down_res@IntersectionSets[["01"]], 
                DEGvenn_down_common=DEGvenn_down_res@IntersectionSets[["11"]], 
                DEGvenn_up_pg.low=DEGvenn_up_res@IntersectionSets[["10"]], 
                DEGvenn_up_pg.high=DEGvenn_up_res@IntersectionSets[["01"]], 
                DEGvenn_up_common=DEGvenn_up_res@IntersectionSets[["11"]]
                )

export(DEG.res,file = '../DEG.res.change.xlsx',row.names=T)
export(pg.pef.tpm,file = '../pg.pef.tpm.xlsx',row.names=T)


## venn for down FC 1.5 for orth
DEGvenn_down = list(pgEpiSC.Low_down = na.omit(pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange < -1.5 ),]$human.genename), 
                    pgEpiSC.High_down = na.omit(pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange < -1.5 ),]$human.genename))
DEGvenn_down_res = Venn(DEGvenn_down)
plot(DEGvenn_down_res,doWeights = F,type="circles")

## venn for up FC 1.5 for orth
DEGvenn_up = list(pgEpiSC.Low_up = na.omit(pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange > 1.5 ),]$human.genename), 
                  pgEpiSC.High_up = na.omit(pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange > 1.5 ),]$human.genename))
DEGvenn_up_res = Venn(DEGvenn_up)
plot(DEGvenn_up_res,doWeights = F,type="circles")

DEG.res <- list(DEGvenn_down_pg.low=DEGvenn_down_res@IntersectionSets[["10"]], 
                DEGvenn_down_pg.high=DEGvenn_down_res@IntersectionSets[["01"]], 
                DEGvenn_down_common=DEGvenn_down_res@IntersectionSets[["11"]], 
                DEGvenn_up_pg.low=DEGvenn_up_res@IntersectionSets[["10"]], 
                DEGvenn_up_pg.high=DEGvenn_up_res@IntersectionSets[["01"]], 
                DEGvenn_up_common=DEGvenn_up_res@IntersectionSets[["11"]]
)

export(DEG.res,file = '../DEG.res.change.orth.xlsx',row.names=T)

pheatmap(t(scale(t(pg.pef.tpm[unique(pg.low_pef.sig$geneid,pg.high_pef.sig$geneid),c(2:9)]))),border_color = NA, clustering_method = 'ward.D',color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),angle_col = '315',show_rownames = F)
pheatmap(
  t(scale(t(pg.pef.tpm[unique(c(pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange < -1.5 ),]$geneid, 
                                pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange < -1.5 ),]$geneid, 
                                pg.low_pef.sig[(pg.low_pef.sig$log2FoldChange > 1.5 ),]$geneid,
                                pg.high_pef.sig[(pg.high_pef.sig$log2FoldChange > 1.5 ),]$geneid)),c(2:9)]))),
           border_color = NA, clustering_method = 'ward.D',color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),angle_col = '315',show_rownames = F)


#----cancer analysis------
library(EnhancedVolcano)
## SF
cancer.SF <- c('SRSF1','ENSSSCG00000036592','SRSF3','SRSF5','SRSF6','SRSF10','TRA2B','ENSSSCG00000000288','ENSSSCG00000036350','HNRNPH1','HNRNPK','HNRNPL','HNRNPM','ENSSSCG00000013421','RBM5','RBM10','RBM39','ESRP1','RBFOX2','QKI')

cancer.SF.tpm <- tpm[cancer.SF,]
cancer.SF.tpm.pg <- pg.pef.tpm[cancer.SF,]

pheatmap(t(scale(t(cancer.SF.tpm.pg[,c(2:9)]))),border_color = NA, clustering_method = 'ward.D',color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=1)(256)),alpha = 1),angle_col = '315',show_rownames = T)


EnhancedVolcano(pg.high_pef,
                lab = pg.high_pef$geneid,
                selectLab = cancer.SF,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1.5,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Cancer SF',
                subtitle = 'P.H vs PEF',
) + coord_flip()


EnhancedVolcano(pg.low_pef,
                lab = pg.low_pef$geneid,
                selectLab = cancer.SF,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1.5,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Cancer SF',
                subtitle = 'P.L vs PEF',
                ) + coord_flip()


## SF downstream


## TF


## downstream genes
can.genes <- c('SERPINE1','	MT2A','UBC','ACTB','TUBA1B','HMGB2','EIF2S3','RPS6','PCNA','STMN1','HMGB1','H2AFZ','VIM','FN1','FOS','DDIT4','NFKBIA','JUNB','TUBB','PKM','LGALS3','BIRC5','TOP2A','NPM1','RPL35','HNRNPH1','SOX4','SOX11','VEGFA','EGFR','MYC','CD44','FOXR2')

tpm.can <- na.omit(pg.pef.tpm[can.genes,-c(1,ncol(pg.pef.tpm))])
pheatmap(t(scale(t(tpm.can))),border_color = NA, clustering_method = 'ward.D2',cluster_cols = F,color = scales::alpha(rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"),alpha=T,bias=.8)(256)),alpha = 1),angle_col = '315',breaks = unique(c(seq(-2,2, length=256))))

EnhancedVolcano(pg.low_pef,
                lab = pg.low_pef$geneid,
                selectLab =can.genes,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1.5,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Cancer genes',
                subtitle = 'P.L vs PEF',
) + coord_flip()

EnhancedVolcano(pg.high_pef,
                lab = pg.high_pef$geneid,
                selectLab =can.genes,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1.5,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Cancer downstream',
                subtitle = 'P.H vs PEF',
) + coord_flip()

## markers from cellmarker
cancer.marker <- c('PROM1','CD44','ABCG2','CD24','CXCR4')


EnhancedVolcano(pg.low_pef,
                lab = pg.low_pef$geneid,
                selectLab = cancer.marker,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1.5,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 2,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'Cancer SF',
                subtitle = 'P.L vs PEF',
) + coord_flip()


EnhancedVolcano(pg.low_pef,
                lab = pg.low_pef$geneid,
                selectLab = plu.gene,
                x = 'log2FoldChange',
                y = 'padj',
                FCcutoff = 1.5,
                pCutoff = 0.05,
                ylab = bquote(~-Log[10]~ 'padj'),
                pointSize = 1.0,
                labSize = 2.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                parseLabels = T,
                col = c('gray80', '#2a9d8f', '#e9c46a', 'red3'),
                colAlpha = 2/5,
                legendPosition = 'bottom',
                legendLabSize = 10,
                legendIconSize = 3.0,
                drawConnectors = TRUE,
                widthConnectors = .5,
                max.overlaps = 300,
                colConnectors = '#14213d',
                title = 'lu genes',
                subtitle = 'P.L vs PEF',
) + coord_flip()

#-----enrichment plot------

## barplot
library(patchwork)
go.bp <- import_list('../pick_GOBP.xlsx')
go.bp[['pg_up']]$Terms <- factor(go.bp[['pg_up']]$Terms, levels = rev(go.bp[['pg_up']]$Terms))
go.bp[['pg_down']]$Terms <- factor(go.bp[['pg_down']]$Terms, levels = rev(go.bp[['pg_down']]$Terms))


scales::show_col(cols)

BLUE <- "#4393C3FF"
RED <- "#D6604DFF"
BLACK <- "#202020"
GREY <- "grey50"


plt <- ggplot(go.bp[['pg_up']]) +
  geom_col(aes(-log10P, Terms), fill = RED, width = 0.6) 
plt

plt <- plt + 
  scale_x_continuous(
    limits = c(0, 20),
    breaks = seq(0, 20, by = 5), 
    expand = c(0, 0), 
    position = "top"  
  ) +
  scale_y_discrete(expand = expansion(add = c(0, 0.5))) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_line(color = "#A8BAC4", size = 0.3),
    axis.ticks.length = unit(0, "mm"),
    axis.title = element_blank(),
    axis.line.y.left = element_line(color = "black"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(family = "Arial", size = 16)
  )

plt

plt <- plt + 
  geom_shadowtext(
    data = go.bp[['pg_up']],
    aes(0, y = Terms, label = Terms),
    hjust = 0,
    nudge_x = 0.3,
    colour = BLACK,
    bg.colour = NA,
    # bg.r = 0.1,
    family = "Arial",
    size = 6
  )

plt.up <- plt

plt

# ## GOplot
# install.packages('GOplot')
# library(GOplot)
# data(EC)
# view(EC$david)
# circ <- circle_dat(EC$david, EC$genelist)
# GOCircle(circ)
# chord <- chord_dat(circ, EC$genes, EC$process)
# reduced_circ <- reduce_overlap(circ, overlap = 0.75)
# 
# GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 5)
# rm(chord,circ,circle_dat,reduced_circ,EC)


## customize

go.bp <- import_list('../GOBP.sort.xlsx')
go.bp.inte <- rbind(go.bp$PEF.GOBP,go.bp$PG.GOBP)
go.bp.inte$Z_score <- (go.bp.inte$Count -mean(go.bp.inte$Count))/sd(go.bp.inte$Count)
go.bp.inte$Type <- c(rep('PEF',20),rep('PG',20))
go.bp.inte$GeneRatio <- go.bp.inte$Count/go.bp.inte$Background
# barplot fail

ggpubr::ggbarplot(go.bp.inte,x = 'Description', y = 'LogQ',
                  fill = 'Type',
                  palette = col[c(1,3)]
                    )


# bubble 
go.bp.inte$Description <- factor(go.bp.inte$Description, levels = rev(go.bp.inte$Description))
ggplot(go.bp.inte[c(1:20),]) + 
  geom_point(aes(GeneRatio,Description,size=Count,color=-LogQ)) +
  scale_color_gradient(low="green",high = "red") +
  labs(color=expression(-log[10](LogQ)),size="Count",  
       x=bquote(~-Log[10]~ 'padj'),y="Terms") +
  theme(legend.title = 'test')

go.bp.inte2 <- go.bp.inte[order(go.bp.inte$GeneRatio),]
go.bp.inte2$Description <- factor(go.bp.inte2$Description,levels = go.bp.inte2$Description)
ggplot(data = go.bp.inte2[1:20,], aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = Count,color = -LogQ)) +
  scale_color_gradient(low = col[1],high = col[3], name=expression(-log[10](padj))) +
  scale_size(range  =  c(0, 6)) +
  labs( y = 'Top 20 GO BP cluster terms') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10)) 

ggplot(data = go.bp.inte2[21:40,], aes(x = GeneRatio, y = Description)) + 
  geom_point(aes(size = Count,color = -LogQ)) +
  scale_fill_brewer(palette = "PuBu", guide = "coloursteps") +
  scale_size(range  =  c(0, 6)) +
  labs( y = 'Top 20 GO BP cluster terms') +
  theme_bw() +
  theme(axis.text.y = element_text(size = 10))



