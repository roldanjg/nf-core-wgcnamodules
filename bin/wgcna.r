#!/usr/bin/env Rscript

library(WGCNA);
library(ComplexHeatmap)
library(circlize)
args = commandArgs(trailingOnly=TRUE)


input_wgcna <- args[1]
allowWGCNAThreads()
message("Input wgcna      (Arg 1): ", input_wgcna)

# Display the current working directory

# Load the WGCNA package

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
femData = read.csv(input_wgcna);


# define the cutting line for the tree of samples
cutthreshold = 1000000 

datExpr0 = as.data.frame(femData[, -1]);
rownames(datExpr0) = femData$gene_id;


gsg = goodSamplesGenes(datExpr0, minNSamples = 2, verbose = 4);
gsg$allOK
if (!gsg$allOK)
{
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
if (sum(!gsg$goodSamples)>0)
printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
# optionally too, uncoment to use, perform a pre-clustering to filter outliers 
# pdf(file = "precut-n-sampleClustering.pdf", width = 36, height = 18);
# sampleTree = hclust(dist(datExpr0), method = "average");
# # Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# # The user should change the dimensions if the window is too large or too small.
# sizeGrWindow(36,18)
# par(cex = 0.6);
# par(mar = c(0,4,2,0))

# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
# cex.axis = 1.5, cex.main = 2)
# # Plot a line to show the cut
# abline(h = cutthreshold, col = "red");
# dev.off()

# clust = cutreeStatic(sampleTree, cutHeight = cutthreshold, minSize = 10)
# table(clust)
# # clust 1 contains the samples we want to keep.
# keepSamples = (clust==1)
# datExpr = datExpr0[keepSamples, ]

datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)


pdf(file = "poscut-n-sampleClustering.pdf", width = 36, height = 18);
sampleTree = hclust(dist(datExpr), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(36,18)
par(cex = 0.6);
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = cutthreshold, col = "red");
dev.off()


powers = c(c(1:10), seq(from = 12, to=30, by=2))
paste(powers)
# https://peterlangfelder.com/2018/11/25/signed-or-unsigned-which-network-type-is-preferable/
# we pick a network type signed:
#  By and large, I recommend using one of the signed varieties, for two main reasons.
#  First, more often than not, direction does matter: it is important to know where 
#  node profiles go up and where they go down, and mixing negatively correlated nodes 
#  together necessarily mixes the two directions together. Second, negatively correlated
#  nodes often belong to different categories
sft = pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = 0.8, verbose = 5, networkType="signed")

pdf(file = "2-n-sft-signed.pdf", width = 9, height = 5);
pdf(file = "nuevo_2-n-sft-signed.pdf", width = 9, height = 5);
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()
paste(sft$powerEstimate)


if (any(is.na(sft$powerEstimate))){
    pwsv = 6
} else {
    pwsv = sft$powerEstimate
}
paste(pwsv)
pwsv = 6
paste(pwsv)

pdf(file = "2-n-softConnectivity-correlation.pdf", width = 9, height = 5);
k=softConnectivity(datE=datExpr,power=pwsv, minNSamples=2)
# Plot a histogram of k and a scale free t
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")
dev.off()
paste(pwsv)
paste("hello")
adjacency = adjacency(datExpr, power = pwsv, type = "signed")
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
sizeGrWindow(12,9)
pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
labels = FALSE, hang = 0.04);

dev.off()

# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
pdf(file = "4-gene_cluster_colors.pdf", width = 12, height = 9);
plotDendroAndColors(geneTree, dynamicMods, "Dynamic Tree Cut",
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05,
main = "Gene dendrogram and module colors")
dev.off()

# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicMods, excludeGrey=TRUE)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
pdf(file = "modules_distance.pdf", width = 12, height = 9);
plot(METree, main = "Clustering of module eigengenes",
xlab = "", sub = "")


MEDissThres = 0.25
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()



# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicMods, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors;
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs;

sizeGrWindow(12, 9)
pdf(file = "modules_dendogram_before_after.pdf", wi = 9, he = 6)
plotDendroAndColors(geneTree, cbind(dynamicMods, mergedColors),
c("Dynamic Tree Cut", "Merged dynamic"),
dendroLabels = FALSE, hang = 0.03,
addGuide = TRUE, guideHang = 0.05)
dev.off()


# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, mergedColors, excludeGrey=TRUE)$eigengenes
MEs = orderMEs(MEs0)

# names (colors) of the modules
modNames = substring(names(MEs), 3)
modNames
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

pdf(file = "ModuleEigenval_distance.biev.pdf", width = 60, height = 60);

mesmats <- data.matrix(MEs)
heatmap(mesmats, margins=c(10,50))
dev.off()

write.csv(mesmats, file = "module_eigengene_values.csv")

# Create the starting data frame
geneInfo0 = data.frame(geneid = names(datExpr),
                      moduleColor = mergedColors)
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, mod], 
                         MMPvalue[, mod]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[mod], sep=""),
                       paste("p.MM.", modNames[mod], sep=""))
}



write.csv(geneInfo0, file = "gene_info.csv")

mesresult = as.matrix(read.csv("module_eigengene_values.csv", header=TRUE, sep=",",row.names = 1,as.is=TRUE))

col_fun1 = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

png("module_trait_relationship.png", res = 300, unit = "cm", height = 15, width = 20)
Heatmap(mesresult, col = col_fun1,
                show_row_names = FALSE,
                column_title = NULL,
                cluster_columns = FALSE,
                cluster_rows = TRUE,
                column_names_rot = 90,
                heatmap_legend_param = list(title ="cor")
                )
      
dev.off()