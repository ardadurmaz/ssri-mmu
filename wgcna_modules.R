library('edgeR')
library('biomaRt')
library('doParallel')
library('WGCNA')
library('fastcluster')

options(stringsAsFactors = FALSE)
expr.count <- read.table('data/pre_processed/batch2_counts_processed.tsv',
                         header = TRUE,
                         sep = '\t',
                         check.names = FALSE)

targets <- data.frame(Gender = c(rep('Female', 18),rep('Male', 18),rep('Embryo', 18)),
                      Tissue = c(rep('Amygdala',3),rep('Hippocampus',3),rep('PrefrontalCortex',3),
                                 rep('Amygdala',3),rep('Hippocampus',3),rep('PrefrontalCortex',3),
                                 rep('Amygdala',3),rep('Hippocampus',3),rep('PrefrontalCortex',3),
                                 rep('Amygdala',3),rep('Hippocampus',3),rep('PrefrontalCortex',3),
                                 rep('Amygdala',3),rep('Hippocampus',3),rep('PrefrontalCortex',3),
                                 rep('Amygdala',3),rep('Hippocampus',3),rep('PrefrontalCortex',3)),
                      Type = c(rep('Control', 9),rep('Treated', 9),
                               rep('Control', 9),rep('Treated', 9),
                               rep('Control', 9),rep('Treated', 9)))
status <- factor(paste(targets$Gender, 
                       targets$Tissue, 
                       targets$Type, 
                       sep = ''))


dge <- DGEList(counts = expr.count, 
               group = status)
dge <- calcNormFactors(dge, 
                       method = 'TMM')
expr <- cpm(dge, 
            normalized.lib.sizes = TRUE, 
            log = TRUE)
expr.norm <- normalizeCyclicLoess(expr)
expr.norm <- t(expr.norm)

enableWGCNAThreads(nThreads = 20)
powers <- seq(from = 1, 
              to = 30, 
              by = 1)
sft <- pickSoftThreshold(expr.norm, 
                         powerVector = powers, 
                         verbose = 5, 
                         networkType = 'unsigned')
plot(sft$fitIndices[,1], -1 * sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     main = paste("Scale independence"), pch = 20)

sft.thr <- 6

adj.data <- adjacency(expr.norm, 
                      type = 'unsigned', 
                      power = sft.thr)
tom.data <- TOMsimilarity(adj.data, 
                          TOMType = 'unsigned')
dissTOM <- 1 - tom.data
saveRDS(dissTOM, file = 'data/TOM-Dissimilarity.rds')

geneTree <- fastcluster::hclust(as.dist(dissTOM), 
                                method = "average")
plot(geneTree, 
     xlab="", 
     sub="", 
     main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, 
                             distM = dissTOM,
                             deepSplit = 4, 
                             pamRespectsDendro = TRUE,
                             minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)

plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")



## Merge Modules ##
MEList <- moduleEigengenes(expr.norm, 
                           colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- fastcluster::hclust(as.dist(MEDiss), 
                              method = 'average')
plot(METree, 
     main = 'Clustering of module eigengenes')
abline(h = 0.25, 
       col = "red")

merge <- mergeCloseModules(expr.norm, 
                          dynamicColors, 
                          cutHeight = 0.25, 
                          verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

plotDendroAndColors(geneTree, 
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleColors = mergedColors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1;

MEs = mergedMEs;
GeneModuleScore <- signedKME(expr.norm, 
                             MEs)

hub.genes <- chooseTopHubInEachModule(datExpr = expr.norm, 
                                      colorh = moduleColors, 
                                      power = 6,
                                      type = 'unsigned')

wgcna.res <- list('MEs' = MEs,
                  'ModuleColors' = moduleColors,
                  'ModuleLabels' = moduleLabels,
                  'GeneModuleScores' = GeneModuleScore,
                  'HubGenes' = hub.genes,
                  'GeneID' = colnames(expr.norm))
saveRDS(wgcna.res,
        file = 'results/WGCNA-Res.rds')


## Generate Plots ##
library(pheatmap)
color.annot <- data.frame('Module' = gsub('^ME', 
                                          replacement = '', 
                                          colnames(wgcna.res$MEs)), 
                          row.names = colnames(wgcna.res$MEs))
pheatmap(cor(wgcna.res$MEs),
         fontsize = 12, 
         scale = 'none',
         main = 'WGCNA Correlation')


