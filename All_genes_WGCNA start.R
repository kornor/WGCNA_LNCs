setwd("~/Bioinformatics Work/Meth & RNA/WGCNA_LNC_ONLY")

### This analysis includes all the genes for mRNA analysis of samples, + lncs
### NOT trimmed to only the methylation related ones

### Load packages
library(WGCNA)
library(flashClust)
library(permute)
library(lattice)
library(vegan)
library(scatterplot3d)
library(MASS)
library(nlme)
library(mgcv)
library(cluster)
library(caret)
library(ape)
library(sparcl)
library(plyr)
library(ggplot2)
library(gplots)
library(RColorBrewer)

## #########################################  This can be skipped - move down to load prepped file
### Load lnc expression file
exp <- read.table( "lnc_for_merge.txt", sep = "\t", header = TRUE, row.names = 1)


## look for genes with missing values
##Exclude genes with results for less than 400 samples; exclude samples with 
#results for less than 20 000 genes. 
gsg = goodSamplesGenes(exp, verbose = 3);
gsg$allOK


## if return is true, all good to go
### Otherwise

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(exp)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(exp)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  exp = exp[gsg$goodSamples, gsg$goodGenes]
}






############ Now need to match and trim the rows (samples) for all sets

### load in the traits and clin

traits <- read.table("Traits_ALL.txt", sep = "\t", header = TRUE, row.names = 1)
clincom <- read.table("Clincom_merged.txt", sep = "\t", header = TRUE, row.names = 1)

list <- intersect(rownames(traits), rownames(clincom))

clincom <- clincom[list,]
exp <- exp[list,]
traits <- traits[list,]



write.table(exp, "Final_exp_lnc.txt", sep = "\t")
write.table(clincom, "Final_clin_lnc.txt", sep = "\t")
write.table(traits, "Final_traits_lnc.txt", sep = "\t")




## Create a dendro of the samples to look for outliers

sampleTree = flashClust(dist(exp), method = "average");

#### Make factor labels for the subtypes
Pam50 <- as.factor(clincom$PAM50Call_RNAseq)
# Relabel blockwise modules


PamColours = labels2colors(Pam50)
count(PamColours)

# Plot the sample tree

sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2, labels = FALSE)

ColorDendrogram(sampleTree, y = PamColours, 
     main = "Sample clustering to detect outliers", 
     sub="", xlab="", 
     labels = FALSE, branchlength = 35)

legend("bottomright",legend=levels(factor(Pam50)),
       fill = (c("turquoise","blue", "brown", "yellow", "green" )), cex = 1)
##########   Different way?

plotDendroAndColors(sampleTree, PamColours,
                    groupLabels = "Pam50 subtype",
                    main = "Sample dendrogram and subtype", dendroLabels = FALSE)
legend("topright",legend=levels(factor(Pam50)),
       fill = (c("turquoise","blue", "brown", "yellow", "green" )), cex = 1)



## Might remove the 4 outliers on the right
## Plot a line to show the cut

abline(h = 1.2, col = "red");

# Determine clusters under the line (use h from above)

#clust = cutreeStatic(sampleTree, cutHeight = 300, minSize = 10)
clust = cutreeStatic(sampleTree, cutHeight = 1.1, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.

keepSamples = (clust==1)
datExpr = exp[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


##Traits

Samples = rownames(datExpr);
traitRows = match(Samples, rownames(traits));
datTraits = traits[traitRows,];


datClin = clincom[traitRows,]

collectGarbage();



###################

# Re-cluster samples to look for underlying patterns

sampleTree2 = flashClust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, 
#grey means missing entry
# need to consider how to do this for discrete variables?

traitColors = numbers2colors(datTraits, 
                             signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap", dendroLabels = FALSE)



#### Note that the side panel - the "basal" group, is hypomethylated across all (comparatively)


### Save out the data and proceed to soft thresholding script ("SingleBlock_moduleCreation.R")


#save(datExpr, datTraits, datClin, file = "lncRNA_total_input.RData")
save(datExpr, datTraits, datClin, file = "lncRNA_trimmed_input.RData")




