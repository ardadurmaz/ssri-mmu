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
wgcna.res <- readRDS('results/WGCNA-Res.rds')

require(ggplot2)
require(reshape2)

for(i in unique(wgcna.res$ModuleColors)){
  local.idx <- which(gsub(pattern = '^kME', replacement = '', colnames(wgcna.res$GeneModuleScores)) == i)
  local.genes <- rownames(wgcna.res$GeneModuleScores)[order(abs(wgcna.res$GeneModuleScores[,local.idx]), decreasing = TRUE)[1:200]]
  local.expr <- expr.norm[,match(local.genes, table = colnames(expr.norm))]
  local.expr.ft <- melt(local.expr)
  colnames(local.expr.ft) <- c('Sample', 'Gene', 'Expression')
  local.expr.ft$Gender <- targets$Gender[match(local.expr.ft$Sample, table = paste('Sample-', 1:54, sep = ''))]
  local.expr.ft$Tissue <- targets$Tissue[match(local.expr.ft$Sample, table = paste('Sample-', 1:54, sep = ''))]
  local.expr.ft$Type <- targets$Type[match(local.expr.ft$Sample, table = paste('Sample-', 1:54, sep = ''))]
  
  p <- ggplot(data = local.expr.ft,
              aes(x = Gender,
                  y = Expression,
                  fill = Type,
                  colour = Type)) +
    geom_point(alpha = 0.2,
               size = 0.5,
               position = position_jitterdodge()) +
    geom_boxplot(alpha = 0.9,
                 notch = TRUE) +
    theme_minimal() +
    ggtitle(sprintf('Module Gene Expression (%s)', i)) +
    theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 16),
          axis.title = element_text(face = 'bold', size = 12),
          axis.text.x = element_text(size = 10),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10))
  ggsave(p, 
         filename = sprintf('plots/wgcna_module_expression_top200__%s.pdf', i),
         width = 10,
         height = 10)
}
