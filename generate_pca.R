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

expr.cpm <- cpm(count.mat, 
                lib.size = colSums(count.mat))
cpm.thr <- (8/median(colSums(count.mat)))*1e6
keep <- rowSums(expr.cpm >= cpm.thr) == ncol(expr.cpm)
count.mat.filt <- count.mat[keep,]

dge <- DGEList(count.mat.filt, group = status)
dge <- calcNormFactors(dge, method = "TMM")
expr.norm <- cpm(dge, normalized.lib.sizes = TRUE, log = TRUE)
expr.norm <- normalizeCyclicLoess(expr.norm)

expr.norm.adult <- expr.norm[,targets$Gender %in% c('Male', 'Female')]
expr.norm.embryo <- expr.norm[,targets$Gender %in% 'Embryo']
targets.adult <- targets[targets$Gender %in% c('Male', 'Female'),]
targets.embryo <- targets[targets$Gender %in% 'Embryo',]

## PCA ##
tissue <- 'Hippocampus'
local.data <- expr.norm.adult[, which(targets.adult$Tissue == tissue)]
local.targets <- targets.adult[targets.adult$Tissue == tissue,]

pca.res <- prcomp(t(local.data), center = TRUE, scale. = TRUE)
plot.data <- data.frame('Coordinate.1' = pca.res$x[,1],
                        'Coordinate.2' = pca.res$x[,2],
                        'Coordinate.3' = pca.res$x[,3],
                        'Gender' = local.targets$Gender,
                        'Tissue' = local.targets$Tissue,
                        'Type' = local.targets$Type)
p <- ggplot(plot.data, 
            aes(x = Coordinate.1, 
                y = Coordinate.2, 
                color = Type,
                shape = Gender)) +
  geom_point(size = 4, alpha = 0.8) +
  ggtitle(sprintf('PCA (Adult-%s)', tissue)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 16))
ggsave(p, filename = sprintf('plots/pca_plot_adult_%s_v2.pdf', tissue), width = 8, height = 8)


## MDS ##
tissue <- 'Hippocampus'
local.data <- expr.norm.adult[, which(targets.adult$Tissue == tissue)]
local.targets <- targets.adult[targets.adult$Tissue == tissue,]

cmd.res <- cmdscale(dist(scale(t(local.data)), method = 'euclidean'))
plot.data <- data.frame('Coordinate.1' = cmd.res[,1],
                        'Coordinate.2' = cmd.res[,2],
                        'Gender' = local.targets$Gender,
                        'Tissue' = local.targets$Tissue,
                        'Type' = local.targets$Type)
p <- ggplot(plot.data, 
            aes(x = Coordinate.1, 
                y = Coordinate.2, 
                color = Type,
                shape = Gender)) +
  geom_point(size = 4, alpha = 0.8) +
  ggtitle(sprintf('MDS (Adult-%s)', tissue)) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = 'bold', size = 16))
ggsave(p, filename = sprintf('plots/mds_plot_adult_%s_v2.pdf', tissue), width = 8, height = 8)
