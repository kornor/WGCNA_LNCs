
###### THIS IS TO COMPARE THE LNC MODS TO THE EXP MODS
### NEED TO DO FOR BOTH TCGA AND META AND THEN PASTE INTO
### THE TRAITS HEATMAPS FILES
################## LET's DO THIS like BRUTUS!

setwd("~/Bioinformatics Work/Meth & RNA/WGCNA_LNC_ONLY")

## Load libraries etc

library(WGCNA)
library(plyr)
library(ggplot2)
library(gplots)

## load in the the lncs files

load(file= "lncRNA_trimmed_input.RData")
load(file = "Dynamic_lncRNA.RData")


### load in the exp file of the modules in all gene set
load( file = "Modules_DS0.RData")
load(file = "MetaAnalysis_trimmed_input.RData")


### ME1A for the match (start with TCGA)

rownames(ME_1A) <- colnames(datExpr1)
  list1 <- intersect(rownames(ME_1A), rownames(datExpr))
    ME <- ME_1A[list1,]
      datExpr <- datExpr[list1,]

other_mod <- ME  # this creates another "traits" file for ExpMods

####### remember that datExpr is now the lncs exp file
### modules1 is the mRNA exp file

nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)


moduleTraitCor = cor(MEs,other_mod, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

### just have a quick look at it
sizeGrWindow(10,6)

par(mar = c(6, 8.5, 3, 3));

# Display the values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(other_mod),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               #textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.85,
               zlim = c(-1,1),
               main = paste("High module expression vs trait genes"))



#### write out this correlation matrix for heatmaps later

write.table(moduleTraitCor, "LNCmods_v_TCGAmods.txt", sep = "\t")
write.table(moduleTraitPvalue, "LNCmods_TCGAmods_pValue.txt", sep = "\t")

#####################################
### Ah yeah crap can't match the META because the samples aren't the same 
### workaround???
