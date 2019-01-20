library(edgeR)
library(ggplot2)

counts <- read.table('data/batch2_counts.tsv',
                     header = TRUE,
                     sep = '\t',
                     comment.char = '#',
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
count.mat <- data.matrix(counts[,7:60])
rownames(count.mat) <- counts$Geneid
colnames(count.mat) <- paste('Sample-', 1:54, sep = '')
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
design <- model.matrix(~0+status)
cm <- makeContrasts(
  FemaleAmygdala = statusFemaleAmygdalaTreated - statusFemaleAmygdalaControl,
  FemaleHippocampus = statusFemaleHippocampusTreated - statusFemaleHippocampusControl,
  FemalePrefrontalCortex = statusFemalePrefrontalCortexTreated - statusFemalePrefrontalCortexControl,
  MaleAmygdala = statusMaleAmygdalaTreated - statusMaleAmygdalaControl,
  MaleHippocampus = statusMaleHippocampusTreated - statusMaleHippocampusControl,
  MalePrefrontalCortex = statusMalePrefrontalCortexTreated - statusMalePrefrontalCortexControl,
  EmbryoAmygdala = statusEmbryoAmygdalaTreated - statusEmbryoAmygdalaControl,
  EmbryoHippocampus = statusEmbryoHippocampusTreated - statusEmbryoHippocampusControl,
  EmbryoPrefrontalCortex = statusEmbryoPrefrontalCortexTreated - statusEmbryoPrefrontalCortexControl,
  AdultAmygdala = ((statusFemaleAmygdalaTreated + statusMaleAmygdalaTreated) / 2) - ((statusFemaleAmygdalaControl + statusMaleAmygdalaControl) / 2),
  AdultHippocampus = ((statusFemaleHippocampusTreated + statusMaleHippocampusTreated) /2) - ((statusFemaleHippocampusControl + statusMaleHippocampusControl)/2),
  AdultPrefrontalCortex = ((statusFemalePrefrontalCortexTreated + statusMalePrefrontalCortexTreated)/2) - ((statusFemalePrefrontalCortexControl + statusMalePrefrontalCortexControl)/2),
  EmbryoComparison = ((statusEmbryoAmygdalaTreated + statusEmbryoHippocampusTreated + statusEmbryoPrefrontalCortexTreated)/3) - ((statusEmbryoAmygdalaControl + statusEmbryoHippocampusControl + statusEmbryoPrefrontalCortexControl)/3),
  AdultComparison = ((statusFemaleAmygdalaTreated + 
                       statusFemaleHippocampusTreated + 
                       statusFemalePrefrontalCortexTreated +
                       statusMaleAmygdalaTreated + 
                       statusMaleHippocampusTreated + 
                       statusMalePrefrontalCortexTreated)/6) - 
    ((statusFemaleAmygdalaControl + 
       statusFemaleHippocampusControl + 
       statusFemalePrefrontalCortexControl +
       statusMaleAmygdalaControl + 
       statusMaleHippocampusControl +
       statusMalePrefrontalCortexControl)/6),
  levels = design)

expr.cpm <- cpm(count.mat, 
                lib.size = colSums(count.mat))
cpm.thr <- (8/median(colSums(count.mat)))*1e6
keep <- rowSums(expr.cpm >= cpm.thr) == ncol(expr.cpm)
count.mat.filt <- count.mat[keep,]

write.table(count.mat.filt, file = 'data/pre_processed/batch2_counts_processed.tsv', 
            sep = '\t', 
            col.names = TRUE, 
            row.names = TRUE, 
            append = FALSE, 
            quote = FALSE)

dge <- DGEList(count.mat.filt, group = status)
dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateGLMRobustDisp(dge, 
                             design, 
                             verbose = TRUE, 
                             maxit = 12)
fit.tagwise <- glmFit(dge, design)


require(biomaRt)
ensembl <- useMart('ENSEMBL_MART_ENSEMBL',
                   dataset = 'mmusculus_gene_ensembl')
mapping <- getBM(mart = ensembl,
                 attributes = c('ensembl_gene_id', 'mgi_symbol'),
                 filters = 'ensembl_gene_id',
                 values = rownames(dge))

deg.list <- sapply(1:ncol(cm), simplify = FALSE, function(c){
  glm.res <- glmLRT(fit.tagwise, contrast = cm[,c])
  results.tagwise <- topTags(glm.res, n = Inf, adjust.method = 'fdr', p.value = 1)$table
  results.tagwise$symbol <- mapping$mgi_symbol[match(rownames(results.tagwise), 
                                                     table = mapping$ensembl_gene_id)]
  results.tagwise$id <- rownames(results.tagwise)
  
  write.table(results.tagwise, 
              file = sprintf("results/gene_stats/batch2_%s.tsv", colnames(cm)[c]),
              col.names = TRUE, 
              row.names = FALSE, 
              sep = "\t", 
              append = FALSE, 
              quote = FALSE)
  
  write.table(subset(results.tagwise, results.tagwise$FDR < 0.1 & abs(results.tagwise$logFC) > log2(1.4)), 
              file = sprintf("results/deg/batch2_deg_%s.tsv", colnames(cm)[c]),
              col.names = TRUE, 
              row.names = FALSE, 
              sep = "\t", 
              append = FALSE, 
              quote = FALSE)
  return(results.tagwise)
})
names(deg.list) <- colnames(cm)


## Exon Level
counts <- read.table('data/batch2_exon_counts.tsv',
                     header = TRUE,
                     sep = '\t',
                     comment.char = '#',
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
idx <- duplicated(paste(counts$Geneid,
                        counts$Start,
                        counts$End, sep = '_'))
counts <- counts[!idx,]


count.mat <- data.matrix(counts[,7:60])
rownames(count.mat) <- paste(counts$Geneid, 
                             paste(counts$Start, counts$End, sep = '-'), 
                             sep = '_')
colnames(count.mat) <- paste('Sample-', 1:54, sep = '')

expr.cpm <- cpm(count.mat, 
                lib.size = colSums(count.mat))
cpm.thr <- (8/median(colSums(count.mat)))*1e6
keep <- rowSums(expr.cpm >= cpm.thr) == ncol(expr.cpm)
count.mat.filt <- count.mat[keep,]
dge <- DGEList(count.mat.filt, group = status)
dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateGLMRobustDisp(dge, design, verbose = TRUE, maxit = 20)
fit.tagwise <- glmFit(dge, design)

library(biomaRt)
ensembl <- useMart('ENSEMBL_MART_ENSEMBL',
                   dataset = 'mmusculus_gene_ensembl')


alt.splice.res <- sapply(1:ncol(cm), simplify = FALSE, function(c){
  sp <- diffSpliceDGE(fit.tagwise, 
                      contrast = cm[,c], 
                      geneid = gsub(pattern = '_.+$', 
                                    replacement = '', 
                                    rownames(dge)),
                      exonid = gsub(pattern = '^ENSMUSG\\d+_',
                                    replacement = '',
                                    gsub(pattern = '-\\d+$',
                                         replacement = '',
                                         rownames(dge))))
  alt.splice <-topSpliceDGE(sp, 
                            test = "Simes", 
                            n = Inf, 
                            FDR = 1)
  mapping <- getBM(mart = ensembl,
                   attributes = c('ensembl_gene_id', 'mgi_symbol'),
                   filters = 'ensembl_gene_id',
                   values = alt.splice$GeneID)
  alt.splice$Symbol <- mapping$mgi_symbol[match(alt.splice$GeneID,
                                                table = mapping$ensembl_gene_id)]
  write.table(subset(alt.splice, alt.splice$FDR < 0.1), 
              file = sprintf("results/alt_splice/batch2_altsplice_%s.tsv", colnames(cm)[c]),
              col.names = TRUE, 
              row.names = FALSE, 
              sep = "\t", 
              append = FALSE, 
              quote = FALSE)
})




plotSpliceDGE(sp, geneid = alt.splice$GeneID[5], genecol = "GeneID")


## GSVA ##
library(GSVA)
library(pheatmap)
library(Rtsne)
library(ggplot2)
reactome <- read.table('data/Ensembl2Reactome_All_Levels.txt',
                       header = FALSE,
                       sep = '\t',
                       quote = '',
                       comment.char = '',
                       stringsAsFactors = FALSE)
reactome <- subset(reactome, reactome$V6 == 'Mus musculus')
reactome.pathways <- split(reactome$V1, reactome$V4)
reactome.pathways <- lapply(reactome.pathways, function(x){
  return(unique(x))
})
reactome.pathways <- reactome.pathways[sapply(reactome.pathways, length) > 8]
expr <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
gsva.scores <- gsva(expr, reactome.pathways)
rownames(targets) <- paste('Sample-', 1:54, sep = '')
pheatmap(gsva.scores, fontsize_row = 1, fontsize_col = 2, annotation_col = targets)


tsne.res <- Rtsne(scale(t(expr)), perplexity = 5, pca_center = FALSE, pca_scale = FALSE)
plot.data <- data.frame('Coordinate.1' = tsne.res$Y[,1],
                        'Coordinate.2' = tsne.res$Y[,2],
                        'Gender' = targets$Gender,
                        'Tissue' = targets$Tissue,
                        'Type' = targets$Type)
ggplot(plot.data, aes(x = Coordinate.1, y = Coordinate.2, color = Type, shape = Gender)) +
  geom_point() +
  theme_minimal() +
  facet_wrap(~Tissue)
