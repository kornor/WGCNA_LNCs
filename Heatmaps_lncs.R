
library(pheatmap)
library(ggplot2)
library(RColorBrewer)


cormat <- read.table("LncMod_Meth_cor_trimmed.txt", sep = "\t", header = TRUE, 
                     row.names = 1)

#### information about par
### mai c(bottom, left, top, right) - the margin in inches
##  mar is as mai but giving the value in lines
### oma c(bottom, left, top, right) of outer margins 
par(oma= c(3,1,1,1))

### Info about heatmap.2
### must be done on a matrix not a dataframe
### dendro can be set to none or some
### key gives you a colour key
### srtCol/ srtRow gives an angle for the text
##sepCol, etc gives a cell outline
### RowSideColors - can set colours for marking the cols or rows



### colorpanel in gplots
###redgreen(n)
###redblue(n)
###bluered(n)
###colorpanel(n, low, mid, high)
## or can set the colour panel before
##eg palette <- colorRampPalette(c('#f0f3ff','#0033BB'))(256)


heatmap.2(as.matrix(t(cormat)),dendrogram="none", 
          notecol="black",col=bluered(200),scale="none",
          key=TRUE, keysize=1.5,key.title = NA, key.xlab = "Correlation",
          density.info="none", trace="none",
          cexRow=1.2,cexCol=1.2, srtCol = 75)



pdf("Heat_test_meth_lnc.pdf",height=8,width=8)
fontsize = 10

pheatmap(t(cormat), color = colorRampPalette(rev(brewer.pal(n = 10, name = "RdYlBu")))(50), 
         kmeans_k = NA, breaks = NA, border_color = "grey60",
         cellwidth = NA, cellheight = NA, scale = "none", 
         cluster_rows = TRUE, cluster_cols = FALSE, fontsize_row = 10,
         annotation_legend = TRUE, drop_levels = TRUE, 
         show_rownames = F,show_colnames = T, 
         treeheight_col = 0,
         main = "Methylation of gene promoters by module")
dev.off()


pdf("Heat_test_meth_lnc1.pdf",height=5,width=3)
par(oma= c(3,1,1,1))
heatmap.2(as.matrix(t(cormat)),dendrogram="row", 
          notecol="black",col=bluered(200),
          scale="none",key=TRUE, keysize=1,
          density.info="none", trace="none",
          labRow = NA, cexRow=1, srtCol = 75, cexCol = 1,
          main = NULL)
dev.off()




tcgaMods <- read.table("ME_tcga.txt", sep = "\t", header = TRUE, row.names = 1)

list <- intersect(row.names(tcgaMods), rownames(datExpr))
tcgaMods <- tcgaMods[list,]
tcgaMods <- tcgaMods[,1:13]
datExpr1 <- datExpr[list,]

nGenes = ncol(datExpr1);
nSamples = nrow(datExpr1);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr1, moduleColors)$eigengenes
MEs = orderMEs(MEs0)



moduleTraitCor = cor(MEs, tcgaMods, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

write.table(moduleTraitCor, "Lnc_vs_TCGA_mods.txt", sep = "\t")
write.table(moduleTraitPvalue, "Lnc_vs_TCGA_mods_pvalue.txt", sep = "\t")


cormat_lnc <- read.table("Lnc_vs_TCGA_mods_cor_trimmed.txt",
                         sep = "\t", header = TRUE,row.names = 1 )

pdf("Heat_test_lncMOD.pdf",height=5,width=8)
par(oma= c(1,1,1,3))

heatmap.2(as.matrix(cormat_lnc),dendrogram="none",
          Rowv = FALSE, Colv = FALSE,
          notecol="black",col=bluered(200),scale="none",
          key=TRUE, keysize=1.5,key.title = NA, key.xlab = "Correlation",
          density.info="none", trace="none",
          cexRow=1.2,cexCol=1.2, srtCol = 75,
          rowsep = 1:nrow(cormat_lnc), colsep = 1:ncol(cormat_lnc), sepcolor = "grey")
dev.off()

##trim out unwanted mods from this

cormat1 <- cormat[1:3,]

pdf("Heat_test_lncMETH.pdf",height=8,width=3)
par(oma= c(1,1,1,3))

heatmap.2(as.matrix(t(cormat1)),dendrogram="row",
          Colv = FALSE, labRow = FALSE,
          notecol="black",col=bluered(200),scale="none",
          key=TRUE, keysize=1.5,key.title = NA, key.xlab = "Correlation",
          density.info="none", trace="none",
          cexRow=1.2,cexCol=1.2, srtCol = 75
)
dev.off()
