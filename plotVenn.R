library(ggplot2)
library(ggforce)

## Read DEG ##
gene.stats <- sapply(list.files(pattern = 'batch2',
                                path = 'results/deg', 
                                full.names = TRUE),
                     simplify = FALSE,
                     function(x){
                       local.data <- read.table(x,
                                                header = TRUE,
                                                sep = '\t',
                                                check.names = FALSE)
                       return(unique(na.omit(local.data$symbol)))
                     })

## Format Names ##
names(gene.stats) <- gsub(pattern = '^batch2_deg_',
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
local.tissue <- gene.stats[grep(pattern = 'hippocampus', ignore.case = TRUE, names(gene.stats))]
match.counts <- matrix(0, ncol = 3, nrow = 8)
res <- apply(as.matrix(combn(0:3, 2)), 2, function(x){
  if(sum(x > 0) == 2){
    length(intersect(local.tissue[[x[1]]],
                     local.tissue[[x[2]]]))  
  }else{
    length(local.tissue[[x[which(x > 0)]]])
  }
})
res <- c(res, length(intersect(local.tissue$EmbryoHippocampus, 
                               intersect(local.tissue$FemaleHippocampus, local.tissue$MaleHippocampus))))


plot.data <- data.frame('Coordinate.1' = c(0, -0.5, 0.5),
                        'Coordinate.2' = c(1, 0, 0),
                        'Group' = c('Embryo', 'Female', 'Male'))
ggplot(plot.data, aes(x0 = Coordinate.1, y0 = Coordinate.2, r = 1, fill = Group)) +
  geom_circle(alpha = 0.6, size = 1, colour = 'grey') +
  coord_fixed() +
  theme_void()
