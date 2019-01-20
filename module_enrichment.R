library(biomaRt)

wgcna.res <- readRDS('results/WGCNA-Res.rds')
module.genes <- split(wgcna.res$GeneID, wgcna.res$ModuleColors, sep = '')
load("data/Curated-GeneList")

## Map ids ##
all.genes <- Reduce(union, curated.genes)
hsa.ensembl <- useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
mapping <- getBM(mart = hsa.ensembl,
                 attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                 filters = 'hgnc_symbol',
                 values = all.genes)
mmu.ensembl <- useMart('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')
orth.mapping <- getLDS(mart = hsa.ensembl,
                       attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                       filters = 'ensembl_gene_id',
                       values = mapping$ensembl_gene_id,
                       attributesL = c('ensembl_gene_id', 'mgi_symbol'),
                       martL = mmu.ensembl)
curated.genes.mapped <- sapply(curated.genes,
                               simplify = FALSE,
                               function(x){
                                 local.ids <- unique(na.omit(mapping$ensembl_gene_id[match(x, table = mapping$hgnc_symbol)]))
                                 local.ids.mmu <- unique(na.omit(orth.mapping$Gene.stable.ID.1[match(local.ids, orth.mapping$Gene.stable.ID)]))
                                 return(local.ids.mmu[local.ids.mmu != ''])
                               })
enr.res <- sapply(module.genes, simplify = TRUE, function(mod){
  local.enr.res <- sapply(curated.genes.mapped, function(cur){
    cur <- cur[cur %in% wgcna.res$GeneID]
    test.mat <- matrix(c(sum(mod %in% cur),
                         sum(!(mod %in% cur)),
                         sum(wgcna.res$GeneID %in% cur),
                         sum(!(wgcna.res$GeneID %in% cur))),
                       byrow = FALSE,
                       ncol = 2)
    local.res <- fisher.test(test.mat, alternative = 'greater')
    return(local.res$p.value)
  })
  return(local.enr.res)
})

require(reshape2)
library(ggplot2)
rownames(enr.res) <- c("SFARI ASD","SFARI ASD (No ID)","SFARI ASD & ID","ID (No SFARI ASD)","ID")
enr.res.ft <- melt(enr.res)
colnames(enr.res.ft) <- c('gene.set', 'module', 'p.value')
enr.res.ft$p.value <- -log10(enr.res.ft$p.value)
enr.res.ft$p.value[enr.res.ft$p.value <= 1] <- NA

p <- ggplot(enr.res.ft, aes(module, gene.set)) + 
  geom_tile(aes(fill = p.value), colour = "white") + 
  scale_fill_gradient2(low = "white", high = "firebrick", midpoint = 0.3, na.value = "white") +
  geom_text(aes(fill = enr.res.ft$p.value, 
                label = round(enr.res.ft$p.value, 2))) +
  theme(panel.background = element_rect(fill = "white"),
        legend.position = "none",
        plot.title = element_text(size = 16, hjust = 0.5),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, vjust = 0.8)) +
  ggtitle("WGCNA Module Enrichment")

temp <- as.data.frame(sapply(module.genes, length))
colnames(temp)[1] <- 'Size'
temp$Module <- rownames(temp)

ggplot(temp, aes(x = Module, y = log(Size))) +
  geom_col(fill = scales::muted('green')) +
  ggtitle('WGCNA Module Size') + 
  ylab('Gene Count (log scale)') + xlab('Module') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 12, face = 'bold'),
        axis.text.y = element_text(size = 12, face = 'bold'),
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 16))
    
