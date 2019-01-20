require(RcisTarget)

gene.stats <- sapply(list.files(pattern = 'batch2',
                                path = 'results/deg', 
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
names(gene.stats) <- gsub(pattern = '^batch2_deg_',
                          replacement = '',
                          gsub(pattern = '\\.tsv$',
                               replacement = '',
                               basename(names(gene.stats))))
wgcna.res <- readRDS('results/WGCNA-Res.rds')
wgcna.genes <- split(wgcna.res$GeneID, wgcna.res$ModuleColors)

## Symbol Mapping ##
all.genes <- Reduce(union, lapply(gene.stats, function(x){names(x)}))
all.genes <- Reduce(union, wgcna.genes)

mmu.ensembl <- useMart('ENSEMBL_MART_ENSEMBL', dataset = 'mmusculus_gene_ensembl')
mapping <- getBM(mart = mmu.ensembl,
                 attributes = c('ensembl_gene_id', 'mgi_symbol'),
                 filters = 'ensembl_gene_id',
                 values = all.genes)
mapping <- na.omit(mapping[mapping$ensembl_gene_id != '' & mapping$mgi_symbol != '',])

## Load TF Data ##
data("motifAnnotations_mgi")
motifRankings <- importRankings('data/mm9-tss-centered-10kb-7species.mc9nr.feather')

tf.enrichment <- sapply(1:length(wgcna.genes), simplify = FALSE, function(x){
  message('Processing Idx: ', x)
  local.genes <- wgcna.genes[[x]]
  geneList <- na.omit(unique(mapping$mgi_symbol[match(local.genes, table = mapping$ensembl_gene_id)]))
  tryCatch({
    # 1. Calculate AUC
    motifs_AUC <- calcAUC(geneList, motifRankings)
    
    # 2. Select significant motifs, add TF annotation & format as table
    motifEnrichmentTable <- addMotifAnnotation(motifs_AUC, 
                                               motifAnnot = motifAnnotations_mgi)
    
    # 3. Identify significant genes for each motif
    # (i.e. genes from the gene set in the top of the ranking)
    # Note: Method 'iCisTarget' instead of 'aprox' is more accurate, but slower
    motifEnrichmentTable_wGenes <- addSignificantGenes(motifEnrichmentTable, 
                                                       geneSets = geneList,
                                                       rankings = motifRankings, 
                                                       nCores = 12,
                                                       method = "iCisTarget")
    return(motifEnrichmentTable_wGenes)
  },
  warning = function(war){
    print('!Warning')
    return(NULL)
  },
  error = function(err){
    print('!Error')
    return(NULL)
  })
})
names(tf.enrichment) <- names(wgcna.genes)
tf.enrichment <- tf.enrichment[!sapply(tf.enrichment, is.null)]

## Get TF Genes ##
require(stringi)
tf.genes <- sapply(tf.enrichment, simplify = FALSE, function(x){
  local.genes <- x$TF_highConf
  local.genes <- local.genes[local.genes != '']
  local.genes <- sapply(local.genes, simplify = FALSE, function(g){
    temp <- unlist(stri_extract_all(g, regex = '(.+?)\\(.+?\\)\\.'))
    temp <- sapply(temp, simplify = FALSE, function(g.g){
      trimws(unlist(strsplit(g.g, split = ';', fixed = TRUE)))
    })
    temp <- Reduce(union, temp)
    temp <- sapply(temp, simplify = FALSE, function(g.g){
      as.character(trimws(gsub(pattern = '\\(.+', replacement = '', g.g)))
    })
    temp <- Reduce(union, temp)
  })
  return(Reduce(union, local.genes))
})

saveRDS(tf.genes, file = 'data/wgcna.tf.genes.rds')

## Plot Network ##
require(Matrix)
require(network)
require(sna)
require(ggplot2)
require(GGally)

ppi <- readRDS('data/stringdb.mmu.symbol.high.conf.rds')
for(i in 1:length(tf.genes)){
  idx <- colnames(ppi) %in% tf.genes[[i]]
  local.ppi <- ppi[idx, idx]
  idx <- colSums(local.ppi) > 0
  local.ppi <- local.ppi[idx, idx]
  graph <- network(as.matrix(local.ppi), directed = FALSE, matrix.type = 'adjacency')
  p <- ggnet2(graph,
              node.size = 6,
              node.color = 'dodgerblue',
              alpha = 0.4,
              edge.size = 0.1,
              label = TRUE,
              label.size = 4) +
    ggtitle(sprintf('TF Network (%s)', names(tf.genes)[i]))
  ggsave(p, filename = sprintf('wgcna_tf_network_%s.pdf', names(tf.genes)[i]), width = 12, height = 12)
}

## Common TF ##
local.common.genes <- Reduce(union, apply(combn(c('turquoise', 'grey60', 'midnightblue', 'salmon'), m = 2), 
                                          2, 
                                          function(x){
  return(intersect(tf.genes[[x[1]]], tf.genes[[x[2]]]))
}))

local.ppi <- matrix(0, 
                    ncol = 4, 
                    nrow = length(local.common.genes), 
                    dimnames = list(local.common.genes, c('turquoise', 'grey60', 'midnightblue', 'salmon')))
for(i in 1:ncol(local.ppi)){
  local.ppi[which(rownames(local.ppi) %in% tf.genes[[which(names(tf.genes) == colnames(local.ppi)[i])]]),i] <- 1
}

require(reshape2)
temp <- melt(local.ppi) 
temp <- subset(temp, temp$value == 1)

graph <- network(temp, 
                 directed = FALSE, 
                 matrix.type = 'edgelist')
graph %v% 'type' <- ifelse(graph %v% 'vertex.names' %in% c('turquoise', 'grey60', 'midnightblue', 'salmon'), 'Module', 'TF')
local.color <- rep('grey', length(graph %v% 'vertex.names'))
local.color[which(graph %v% 'type' == 'Module')] <- get.vertex.attribute(graph, 'vertex.names')[which(graph %v% 'type' == 'Module')]
graph %v% 'color' <- local.color



p <- ggnet2(graph,
            node.size = 6,
            node.color = local.color,
            alpha = 0.8,
            edge.size = 0.1,
            label = TRUE,
            label.size = 4) +
  ggtitle(sprintf('TF Network', names(tf.genes)[i]))
