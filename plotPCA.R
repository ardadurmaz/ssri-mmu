library('edgeR')
library('biomaRt')
library('fastcluster')
library('ggplot2')

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
expr.norm <- normalizeCyclicLoess(expr,
                                  method = 'pairs')
expr.norm <- scale(t(expr.norm))

pca.res <- prcomp(expr.norm,
                  center = FALSE,
                  scale. = FALSE)

plot.data <- data.frame('Coordinate.1' = pca.res$x[,1],
                        'Coordinate.2' = pca.res$x[,2],
                        'Coordinate.3' = pca.res$x[,3],
                        'Coordinate.4' = pca.res$x[,4],
                        'Gender' = targets$Gender,
                        'Tissue' = targets$Tissue,
                        'Type' = targets$Type)
ggplot(plot.data,
       aes(x = Coordinate.1,
           y = Coordinate.2,
           color = Type,
           shape = Gender)) +
  geom_point(size = 4) +
  theme_minimal() + 
  facet_wrap(~Tissue) +
  ggtitle('PCA') +
  theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
        strip.text = element_text(size = 10, face = 'bold'))
