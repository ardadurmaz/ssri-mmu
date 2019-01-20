library('fgsea')

## Read DEG Data ##
gene.stats <- sapply(list.files(pattern = 'batch1',
                                path = 'results/gene_stats', 
                                full.names = TRUE),
                     simplify = FALSE,
                     function(x){
                       local.data <- read.table(x,
                                                header = TRUE,
                                                sep = '\t',
                                                check.names = FALSE)
                       local.stats <- setNames(local.data$logFC, nm = local.data$id)
                       local.stats <- local.stats[order(abs(local.stats), decreasing = TRUE)]
                       local.stats <- local.stats[!duplicated(names(local.stats))]
                       return(sort(local.stats, decreasing = TRUE))
                     })

## Format Names ##
names(gene.stats) <- gsub(pattern = '^batch1_',
                          replacement = '',
                          gsub(pattern = '\\.tsv$',
                               replacement = '',
                               basename(names(gene.stats))))

## Read Pathway Data ##
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
gc()

gsea.res <- sapply(gene.stats, simplify = FALSE, function(x){
  local.res <- fgsea(pathways = reactome.pathways,
                     stats = x, 
                     nperm = 100000, 
                     minSize = 8, 
                     nproc = 20)
  return(local.res)
})

gsea.res.ft <- sapply(gsea.res, simplify = FALSE, function(x){
  x <- as.data.frame(x)
  res <- data.frame('pathway' = x$pathway,
                    'adj.p.val' = x$padj,
                    'nes' = x$NES,
                    'leading.edge' = sapply(x$leadingEdge, function(t){return(paste(t,collapse=','))}))
  return(res)
})

## Plot Enrichment ##
require(ggplot2)
require(reshape2)

gsea.res.ft <- sapply(gsea.res, function(x){
  x <- as.data.frame(x)
  x <- x[match(names(reactome.pathways),
               table = x$pathway),]
  x$pathway <- names(reactome.pathways)
  x$NES[x$padj >= 0.1] <- NA
  return(setNames(x$NES, x$pathway))
})
gsea.res.ft <- gsea.res.ft[!(rowSums(is.na(gsea.res.ft)) == ncol(gsea.res.ft)),]
gsea.res.ft <- melt(gsea.res.ft)
colnames(gsea.res.ft) <- c('Pathway', 'Comparison', 'NES')
gsea.res.ft <- na.omit(gsea.res.ft)

ggplot(gsea.res.ft, aes(x = Comparison, y = Pathway, color = NES)) +
  geom_point() +
  theme_classic() +
  ggtitle('Reactome Pathway Enrichment') +
  theme(axis.text.x = element_text(angle = 60, 
                                   hjust = 1,
                                   size = 12),
        axis.text.y = element_text(size = 3),
        axis.ticks = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 16),
        panel.border = element_blank()) +
  scale_color_gradient2(low = 'dodgerblue',
                        mid = 'white',
                        high = 'firebrick')

## Module Enrichment ##
wgcna.res <- readRDS('~/MusMusculus-RNASeq/results/WGCNA-Res.rds')
wgcna.modules <- split(wgcna.res$GeneID, wgcna.res$ModuleColors)
all.genes <- Reduce(union, wgcna.modules)

module.enr <- sapply(wgcna.modules, function(m){
  local.enr <- sapply(reactome.pathways, function(p){
    test.mat <- matrix(c(sum(m %in% p),
                         sum(!(m %in% p)),
                         sum(all.genes %in% p),
                         sum(!(all.genes %in% p))),
                       byrow = FALSE, 
                       ncol = 2)
    test.res <- fisher.test(test.mat, alternative = 'greater')
    return(test.res$p.value)
  })
  local.enr <- p.adjust(local.enr, method = 'bonferroni')
  return(local.enr)
})
module.enr <- module.enr[!apply(module.enr, 1, function(x){all(x >= 0.1)}),]
module.enr.ft <- melt(module.enr)
colnames(module.enr.ft) <- c('Pathway', 'Module', 'P.val')

ggplot(module.enr.ft, 
       aes(x = Module, 
           y = Pathway, 
           color = -log10(P.val))) +
  geom_point() +
  theme_classic() +
  ggtitle('Reactome Pathway Enrichment') +
  theme(axis.text.x = element_text(angle = 60, 
                                   hjust = 1,
                                   size = 12),
        axis.text.y = element_text(size = 8),
        axis.ticks = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 16),
        panel.border = element_blank()) +
  scale_color_gradient2(low = 'white',
                        high = 'firebrick')

## Plot Expression ##
library(edgeR)
library(ggplot2)
library(biomaRt)
library(pheatmap)

count.mat <- read.table('/media/drmz/SSDA/MusMusculus-RNASeq/data/pre_processed/batch1_counts_processed.tsv',
                        header = TRUE,
                        sep = '\t',
                        stringsAsFactors = FALSE,
                        check.names = FALSE)
targets <- c('Control', 'Control', 'Control', 'Treated', 'Treated') 
status <- factor(targets)

dge <- DGEList(count.mat, 
               group = status)
dge <- calcNormFactors(dge, method = "TMM")
expr.norm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)

leading.edge.genes <- unlist(strsplit(as.character(gsea.res.ft$gene_stats$leading.edge[which(gsea.res.ft$gene_stats$pathway == 'Rho GTPase cycle')]), split = ','))
expr.norm.filt <- expr.norm[rownames(expr.norm) %in% leading.edge.genes,]

ensembl <- useMart('ENSEMBL_MART_ENSEMBL',
                   dataset = 'mmusculus_gene_ensembl')
mapping <- getBM(mart = ensembl,
                 attributes = c('ensembl_gene_id', 'mgi_symbol'),
                 filters = 'ensembl_gene_id',
                 values = rownames(expr.norm.filt))
rownames(expr.norm.filt) <- mapping$mgi_symbol[match(rownames(expr.norm.filt), table = mapping$ensembl_gene_id)]
annot <- data.frame('Group' = as.character(status),
                    row.names = colnames(expr.norm.filt))
pheatmap(expr.norm.filt, show_colnames = FALSE, annotation_col = annot, cluster_cols = FALSE, main = 'Rho GTPase Cycle')


## WGCNA Enrichment ##
wgcna.res <- readRDS('results/WGCNA-Res.rds')
module.genes <- split(wgcna.res$GeneID, wgcna.res$ModuleColors, sep = '')
gene.stats <- sapply(list.files(pattern = 'batch2',
                                path = 'results/gene_stats', 
                                full.names = TRUE),
                     simplify = FALSE,
                     function(x){
                       local.data <- read.table(x,
                                                header = TRUE,
                                                sep = '\t',
                                                check.names = FALSE)
                       local.stats <- setNames(local.data$logFC, nm = local.data$id)
                       local.stats <- local.stats[order(abs(local.stats), decreasing = TRUE)]
                       local.stats <- local.stats[!duplicated(names(local.stats))]
                       return(sort(local.stats, decreasing = TRUE))
                     })

## Format Names ##
names(gene.stats) <- gsub(pattern = '^batch2_',
                          replacement = '',
                          gsub(pattern = '\\.tsv$',
                               replacement = '',
                               basename(names(gene.stats))))
gene.stats <- gene.stats[!(names(gene.stats) %in% c('AdultAmygdala', 
                                                    'AdultComparison',
                                                    'AdultHippocampus',
                                                    'AdultPrefrontalCortex',
                                                    'Embryo',
                                                    'EmbryoComparison'))]
gsea.res <- sapply(gene.stats, simplify = FALSE, function(x){
  local.res <- fgsea(pathways = module.genes,
                     stats = x, 
                     nperm = 100000, 
                     minSize = 8, 
                     nproc = 20)
  return(local.res)
})

gsea.res.ft <- sapply(gsea.res, simplify = FALSE, function(x){
  x <- as.data.frame(x)
  res <- data.frame('pathway' = x$pathway,
                    'adj.p.val' = x$padj,
                    'nes' = x$NES,
                    'leading.edge' = sapply(x$leadingEdge, function(t){return(paste(t,collapse=','))}))
  return(res)
})

all.paths <- as.character(gsea.res.ft$EmbryoAmygdala$pathway)
gsea.res.ft.ft <- sapply(gsea.res.ft, function(x){
  local <- setNames(x$nes, nm = x$pathway)
  local[x$adj.p.val >= 0.1] <- NA
  local <- local[match(all.paths, table = names(local))]
  return(local)
})
gsea.res.ft.ft <- melt(gsea.res.ft.ft)
colnames(gsea.res.ft.ft) <- c('Module', 'Comparison', 'NES')

ggplot(gsea.res.ft.ft, 
       aes(x = Module, 
           y = Comparison)) +
  geom_tile(aes(fill = NES)) +
  ggtitle('WGCNA Module Enrichment') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 16)) +
  scale_fill_gradient2(low = scales::muted('dodgerblue'), 
                        high = scales::muted('firebrick'), 
                        na.value = 'white')
##
sig.modules <- c('magenta', 'brown')
temp <- subset(gsea.res.ft.ft, gsea.res.ft.ft$Module == 'magenta')
temp$Gender <- rep("", nrow(temp))
temp$Gender[grep(pattern = 'Male', temp$Comparison)] <- 'Male'
temp$Gender[grep(pattern = 'Female', temp$Comparison)] <- 'Female'
temp$Gender[grep(pattern = 'Embryo', temp$Comparison)] <- 'Embryo'

temp$Tissue <- rep("", nrow(temp))
temp$Tissue[grep(pattern = 'Amygdala', temp$Comparison)] <- 'Amygdala'
temp$Tissue[grep(pattern = 'Hippocampus', temp$Comparison)] <- 'Hippocampus'
temp$Tissue[grep(pattern = 'PrefrontalCortex', temp$Comparison)] <- 'PrefrontalCortex'


ggplot(temp, 
       aes(x = Tissue, 
           y = Gender)) +
  geom_tile(aes(fill = NES)) +
  geom_text(aes(label = round(NES, digits = 2))) +
  ggtitle('WGCNA Module Enrichment (magenta)') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.ticks = element_blank(),
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 16)) +
  scale_fill_gradient2(low = scales::muted('dodgerblue'), 
                       high = scales::muted('firebrick'), 
                       na.value = 'white')
