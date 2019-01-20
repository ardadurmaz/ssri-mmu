library(Matrix)
library(biomaRt)

adj <- readRDS('~/Research/resources/stringdb/stringdb.mmu.symbol.med.conf.rds')
wgcna.res <- readRDS('/media/drmz/SSDA/MusMusculus-RNASeq/results/WGCNA-Res.rds')
module.genes <- split(wgcna.res$GeneID, wgcna.res$ModuleColors, sep = '')
sig.modules <- module.genes[which(names(module.genes)  %in% 'purple')]

mmu.mart <- useMart('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')
mmu.mapping <- getBM(mart = mmu.mart,
                     attributes = c('ensembl_gene_id', 'mgi_symbol'),
                     filters = c('ensembl_gene_id'),
                     values = Reduce(union, sig.modules))

sig.modules.mapped <- lapply(sig.modules, function(x){
  x.mapped <- na.omit(mmu.mapping$mgi_symbol[match(x, table = mmu.mapping$ensembl_gene_id)])
  return(unique(x.mapped[!(x.mapped == '')]))
})

idx <- colnames(adj) %in% Reduce(union, sig.modules.mapped)
local.adj <- adj[idx, idx]

## Plot Network ##
require(network)
require(sna)
require(ggplot2)
require(ggnet)

graph <- network(as.matrix(local.adj), directed = FALSE, matrix.type = 'adjacency')
graph %v% 'Module' <- ifelse(graph %v% 'vertex.names' %in% sig.modules.mapped$purple, 'purple')
ggnet2(graph, 
       color = 'Module',
       label = TRUE,
       size = 6,
       label.size = 1)
