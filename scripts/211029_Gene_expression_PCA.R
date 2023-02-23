library(PCAtools)
library(dplyr)
library(tidyverse)
library(stringr)
library(envalysis)

## Load upstream processed data
load("data/rnaseq/processed/Processed_RNAseq.Rdata")

meta_pca <- meta %>% select(Age, Genotype, Replicate, Group)
rownames(meta_pca) <- paste0(meta_pca$Age, "_", meta_pca$Genotype, "_", meta_pca$Replicate)

chr_list <- read.csv('data/chromosomes/chr_list', colClasses = c("character"))

## Prep data for PCA
pca_mat <- data.frame(cbind(gene.dge.filt$genes, gene.dge.filt$lcpm)) %>% filter(chromosome_name %in% chr_list$X1) %>% 
  select(-chromosome_name, -strand, -start_position, -end_position, -description, -transcript_count, -gene_biotype, -source) 
rownames(pca_mat) <- make.names(pca_mat$external_gene_name, unique = T)
pca_mat <- pca_mat %>% select(-external_gene_name)

## Make PCA
p <- pca(mat = pca_mat, center=TRUE, scale=TRUE, metadata = meta_pca, removeVar = 0.1)

# ## Check variation exxplained per PC
# screeplot(p, hline=50)
# ggsave(filename = paste0("../raw_plots/",format(Sys.time(), "%Y_%m_%d"),"_Gene_expression_scree_plot.pdf"), 
#        device = 'pdf', 
#        height = 7, width = 7, units = "in", dpi = "retina")
# 
# ## Check correlation betwen PCs and factors
# eigencorplot(p, metavars = c("Genotype", "Age", "Replicate"))
# ggsave(filename = paste0("../raw_plots/",format(Sys.time(), "%Y_%m_%d"),"_Gene_expression_eigencor_plot.pdf"), 
#        device = 'pdf', 
#        height = 7, width = 7, units = "in", dpi = "retina")
# 
# 
# ggsave(filename = paste0("../raw_plots/",format(Sys.time(), "%Y_%m_%d"),"_Gene_expression_PCA_plot.pdf"), 
#        device = 'pdf', 
#        height = 7, width = 7, units = "in", dpi = "retina")



## Scree Plot
pscree <- screeplot(p, components = getComponents(p, 1:10),
                    hline = 75, axisLabSize = 14, titleLabSize = 20,
                    returnPlot = FALSE, title='') +
  geom_label(aes(7, 58, label = '75% explained variation', vjust = -1, size = 8))

## Eigencor Plot
peigencor <- eigencorplot(p,
                          components = getComponents(p, 1:5),
                          metavars = c("Genotype", "Age", "Replicate"),
                          cexCorval = 1.0,
                          fontCorval = 2,
                          rotLabX = 45,
                          col = jet
)

## Biplot
pbiplot <- biplot(p, 
                  colby = 'Group', 
                  colkey = c('P28_WT' = "#61ADF1", 'P28_KO' = "#FFE186",'M12_WT' =  "#10518A",'M12_KO' = "#F7B906"),
                  encircle = T)

## Loadings Plot
PC1_loadings <- p$loadings %>% select(PC1) %>% mutate(abs_PC = abs(PC1)) %>% arrange(-abs_PC) %>% top_n(25, abs_PC)
ggplot(data = PC1_loadings, aes(x = "PC1", y = abs_PC)) + geom_point() + geom_text(label=rownames(PC1_loadings))

plotloadings(p, components = getComponents(p, 1:1), rangeRetain = 0.003, absolute = T, labSize = 5)
options(ggrepel.max.overlaps = Inf)
ploadings <- plotloadings(p, rangeRetain = 0.0021, labSize = 4,
                          title = 'Loadings plot', axisLabSize = 12,
                          components = getComponents(p, 1:2),
                          caption = 'Top 10 variables',
                          shape = 24, shapeSizeRange = c(4, 8),
                          col = c('limegreen', 'black', 'red3'),
                          legendPosition = 'none',
                          drawConnectors = T,
                          returnPlot = FALSE)

ppairs <- pairsplot(p, components = getComponents(p, c(1:3)),
                    triangle = TRUE, trianglelabSize = 15,
                    hline = 0, vline = 0,
                    pointSize = 3, gridlines.major = FALSE, gridlines.minor = FALSE,
                    colby = 'Group', 
                    colkey = c('P28_WT' = "#61ADF1", 'P28_KO' = "#FFE186", 'M12_WT' =  "#10518A", 'M12_KO' = "#F7B906"),
                    title = '', plotaxes = FALSE,
                    margingaps = unit(c(0.01, 0.01, 0.01, 0.01), 'cm'),
                    returnPlot = FALSE)

pbiplot <- biplot(p,
                  # loadings parameters
                  showLoadings = TRUE,
                  lengthLoadingsArrowsFactor = 1.5,
                  sizeLoadingsNames = 4,
                  colLoadingsNames = 'red4',
                  # other parameters
                  lab = NULL,
                  colby = 'ER', colkey = c('ER+'='royalblue', 'ER-'='red3'),
                  hline = 0, vline = c(-25, 0, 25),
                  vlineType = c('dotdash', 'solid', 'dashed'),
                  gridlines.major = FALSE, gridlines.minor = FALSE,
                  pointSize = 5,
                  legendPosition = 'none', legendLabSize = 16, legendIconSize = 8.0,
                  shape = 'Grade', shapekey = c('Grade 1'=15, 'Grade 2'=17, 'Grade 3'=8),
                  drawConnectors = FALSE,
                  title = 'PCA bi-plot',
                  subtitle = 'PC1 versus PC2',
                  caption = '27 PCs ≈ 80%',
                  returnPlot = FALSE)





library(cowplot)
library(ggplotify)

top_row <- plot_grid(pscree, ppairs, pbiplot,
                     ncol = 3,
                     labels = c('A', 'B', 'C'),
                     label_fontfamily = 'serif',
                     label_fontface = 'bold',
                     label_size = 22,
                     align = 'h',
                     rel_widths = c(1.10, 0.80, 1.10))

bottom_row <- plot_grid(ploadings,
                        as.grob(peigencor),
                        ncol = 2,
                        labels = c('D', 'E'),
                        label_fontfamily = 'serif',
                        label_fontface = 'bold',
                        label_size = 22,
                        align = 'h',
                        rel_widths = c(0.8, 1.2))

plot_grid(top_row, bottom_row, ncol = 1,
          rel_heights = c(1.1, 0.9))
