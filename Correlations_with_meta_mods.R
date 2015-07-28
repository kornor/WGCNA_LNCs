
### load in the exp file of the modules in all gene set
load( file = "Modules_DS0.RData")
load(file = "MetaAnalysis_trimmed_input.RData")
### ME1A for the match

rownames(ME_1A) <- colnames(datExpr1)

list1 <- intersect(rownames(ME_1A), rownames(datExpr))

ME <- ME_1A[list1,]
datExpr <- datExpr[list1,]

other_mod <- ME


nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)


moduleTraitCor = cor(MEs,other_mod, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");

##########Will only display correlations
textMatrix = paste(signif(moduleTraitCor, 2),sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

# Display the values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(other_mod),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.85,
               zlim = c(-1,1),
               main = paste("High module expression vs trait genes"))



#### write out this correlation matrix for heatmaps later
