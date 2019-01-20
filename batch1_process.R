library(edgeR)
library(ggplot2)
library(biomaRt)

counts <- read.table('data/batch1_counts.tsv',
                     header = TRUE,
                     sep = '\t',
                     comment.char = '#',
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
count.mat <- data.matrix(counts[,7:12])
rownames(count.mat) <- counts$Geneid
colnames(count.mat) <- paste('Sample-', 1:6, sep = '')
rm(counts)

count.mat <- count.mat[,-5]

targets <- c('Control', 'Control', 'Control', 
             'Treated', 'Treated') #,'Treated')
status <- factor(targets)
design <- model.matrix(~0+status)
cm <- makeContrasts(
  TC = statusTreated - statusControl,
  levels = design)

## Filter Counts ##
expr.cpm <- cpm(count.mat, lib.size = colSums(count.mat), log = FALSE)
cpm.thr <- (8/median(colSums(count.mat)))*1e6
keep <- rowSums(expr.cpm >= cpm.thr) == ncol(expr.cpm)
count.mat.filt <- count.mat[keep,]

## Filter Protein Coding
# rownames(mapping) <- 1:nrow(mapping)
# idx <- mapping$gene_biotype == 'protein_coding'
# count.mat.filt <- count.mat.filt[idx,]
write.table(count.mat.filt, 
            file = 'data/pre_processed/batch1_counts_processed.tsv',
            col.names = TRUE,
            row.names = TRUE,
            sep = '\t',
            quote = FALSE,
            append = FALSE)

dge <- DGEList(count.mat.filt, group = status)
dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateGLMRobustDisp(dge, 
                             design = design, 
                             maxit = 12, 
                             verbose = TRUE)
fit <- glmFit(dge, design = design)
glm.res <- glmLRT(fit, contrast = cm[, 1])
results.tagwise <- topTags(glm.res, 
                           n = Inf, 
                           adjust.method = 'fdr', 
                           p.value = 1)$table

ensembl <- useMart('ENSEMBL_MART_ENSEMBL',
                   dataset = 'mmusculus_gene_ensembl')
mapping <- getBM(mart = ensembl,
                 attributes = c('ensembl_gene_id', 'mgi_symbol'),
                 filters = 'ensembl_gene_id',
                 values = rownames(results.tagwise))
results.tagwise$symbol <- mapping$mgi_symbol[match(rownames(results.tagwise),
                                                   table = mapping$ensembl_gene_id)]
results.tagwise$id <- rownames(results.tagwise)
results.tagwise.sig <- subset(results.tagwise, results.tagwise$FDR < 0.1 & abs(results.tagwise$logFC) > log2(1.4))

write.table(results.tagwise.sig,
            file = 'results/deg/batch1_deg.tsv',
            col.names = TRUE,
            row.names = FALSE,
            sep = '\t',
            append = FALSE,
            quote = FALSE)

write.table(results.tagwise,
            file = 'results/gene_stats/batch1_gene_stats.tsv',
            col.names = TRUE,
            row.names = FALSE,
            sep = '\t',
            append = FALSE,
            quote = FALSE)



## Exon Level
counts <- read.table('data/batch1_exon_counts.tsv',
                     header = TRUE,
                     sep = '\t',
                     comment.char = '#',
                     check.names = FALSE,
                     stringsAsFactors = FALSE)
idx <- duplicated(paste(counts$Geneid,
                        counts$Start,
                        counts$End, sep = '_'))
counts <- counts[!idx,]
count.mat <- data.matrix(counts[,7:12])
rownames(count.mat) <- paste(counts$Geneid, 
                             paste(counts$Start, 
                                   counts$End, 
                                   sep = '-'), 
                             sep = '_')
colnames(count.mat) <- paste('Sample-', 1:6, sep = '')
targets <- c('Control', 'Control', 'Control', 'Treated', 'Treated', 'Treated')
status <- factor(targets)
design <- model.matrix(~0+status)
cm <- makeContrasts(
  TC = statusTreated - statusControl,
  levels = design)

expr.cpm <- cpm(count.mat, lib.size = colSums(count.mat))
cpm.thr <- (12/median(colSums(count.mat)))*1e6
keep <- rowSums(expr.cpm >= cpm.thr) == ncol(expr.cpm)
count.mat.filt <- count.mat[keep,]
dge <- DGEList(count.mat.filt, group = status)
dge <- calcNormFactors(dge, method = "TMM")
dge <- estimateGLMRobustDisp(dge, 
                             design = design, 
                             maxit = 20, 
                             verbose = TRUE)
fit <- glmFit(dge, design = design)
sp <- diffSpliceDGE(fit, 
                    contrast = cm[,1], 
                    geneid = gsub(pattern = '_.+$', 
                                  replacement = '', 
                                  rownames(dge)),
                    exonid = gsub(pattern = '^ENSMUSG\\d+_',
                                  replacement = '',
                                  rownames(dge)))

alt.splice <-topSpliceDGE(sp, 
                          test = "Simes", 
                          n = Inf, 
                          FDR = 0.1)
ensembl <- useMart('ENSEMBL_MART_ENSEMBL',
                   dataset = 'mmusculus_gene_ensembl')
mapping <- getBM(mart = ensembl,
                 attributes = c('ensembl_gene_id', 'mgi_symbol'),
                 filters = 'ensembl_gene_id',
                 values = alt.splice$GeneID)
alt.splice$Symbol <- mapping$mgi_symbol[match(alt.splice$GeneID,
                                              table = mapping$ensembl_gene_id)]
write.table(alt.splice,
            file = 'results/alt_splice/batch1_alt_splice.tsv',
            col.names = TRUE,
            row.names = FALSE,
            append = FALSE,
            quote = FALSE,
            sep = '\t')
