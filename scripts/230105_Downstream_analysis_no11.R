## 01-05-2023
## Downstream analysis of gene expression changes in PILRB1 KO mice
## And figure generation

## Load the upstream processed data
load("data/rnaseq/processed_no11/Processed_RNAseq_w_CPM.Rdata")

## load libraries
library(gprofiler2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)
library(dplyr)

## Set focal comparison 
df <- DE$P28_WTvKO  ## P28 mice WT vs KO 
df <- DE$M12_WTvKO  ## 12M mice WT vs KO
#df <- DE$P28_v_M12  ## P28 vs 12M (regardless of genotype)
#df <- DE$WT_v_KO    ## WT vs KO (regardless of age)

## Get log2 fold change 
gene_list <- df$logFC

## name the vector
names(gene_list) <- df$ensembl_ID

## Sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


## Run GSEA on GO terms
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")


## Up Only
gse_up <- clusterProfiler::filter(gse, NES>0)
gse_up <- setReadable(gse_up, OrgDb = org.Mm.eg.db)
gse_up_sim <- pairwise_termsim(gse_up)
p <- treeplot(gse_up_sim,  hclust_method = "average")
p

p$data %>% filter(isTip == T) %>% group_by(group) %>% slice_head(n=1)

selected_pathways <- p$data %>% filter(group %in% c("cluster_1", "cluster_2", "cluster_4")) %>% filter(!is.na(label)) %>% pull(label) ## Olfactory
clusterProfiler::heatplot(gse_up, showCategory = selected_pathways, foldChange = gene_list) + scale_x_discrete(guide = guide_axis(n.dodge = 2))

selected_pathways <- p$data %>% filter(group %in% c("cluster_5")) %>% filter(!is.na(label)) %>% pull(label) ## Pheromone
clusterProfiler::heatplot(gse_up, showCategory = selected_pathways, foldChange = gene_list) + scale_x_discrete(guide = guide_axis(n.dodge = 1))

selected_pathways <- p$data %>% filter(group %in% c("cluster_3")) %>% filter(!is.na(label)) %>% pull(label) ## G prot
clusterProfiler::heatplot(gse_up, showCategory = selected_pathways, foldChange = gene_list) + scale_x_discrete(guide = guide_axis(n.dodge = 2))




## Down Only
gse_down <- clusterProfiler::filter(gse, NES<0)
gse_down <- setReadable(gse_down, OrgDb = org.Mm.eg.db)
gse_down_sim <- pairwise_termsim(gse_down)
p <- treeplot(gse_down_sim, hclust_method = "average")
p

p$data %>% filter(isTip == T) %>% group_by(group) %>% slice_head(n=1)

selected_pathways <- p$data %>% filter(group %in% c("cluster_2", "cluster_4")) %>% filter(!is.na(label)) %>% pull(label) ## Ribosome
clusterProfiler::heatplot(gse_down, showCategory = selected_pathways, foldChange = gene_list) + scale_x_discrete(guide = guide_axis(n.dodge = 4))

selected_pathways <- p$data %>% filter(group %in% c("cluster_3", "cluster_5")) %>% filter(!is.na(label)) %>% pull(label) ## Mitochondrial
clusterProfiler::heatplot(gse_down, showCategory = selected_pathways, foldChange = gene_list) + scale_x_discrete(guide = guide_axis(n.dodge = 4))

selected_pathways <- p$data %>% filter(group %in% c("cluster_1")) %>% filter(!is.na(label)) %>% pull(label) ## Lens
clusterProfiler::heatplot(gse_down, showCategory = selected_pathways, foldChange = gene_list) + scale_x_discrete(guide = guide_axis(n.dodge = 1))


## Convert gene IDs for gseKEGG function
## We will lose some genes here because not all IDs will be converted
ids <- bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)

## Remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids <- ids[!duplicated(ids[c("ENSEMBL")]),]

## Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 <- df[df$ensembl_ID %in% dedup_ids$ENSEMBL,]

## Create a new column in df2 with the corresponding ENTREZ IDs
df2$ENTREZID <- dedup_ids$ENTREZID

## Create a vector of the gene universe
kegg_gene_list <- df2$logFC

## Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$ENTREZID

## Sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

##  Run GSEA on KEGG terms
kegg_organism = "mmu"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")

## KEGG term Dotplot
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

library("pathview")
mmu00590 <- pathview(gene.data  = kegg_gene_list,
                     pathway.id = "mmu00590",
                     species    = "mmu",
                     limit      = list(gene=max(abs(kegg_gene_list)), cpd=1),
                    low = "blue", high = "red")
knitr::include_graphics("mmu04740.pathview.png")


## KEGG term GSEA plot
## Look at the list of GSEA gene sets and pick one to plot by ID
## You can check all available sets by uncommenting the next line
# gse$Description 
geneset = 1
gseaplot(kk2, by = "all", title = kk2$Description[geneset], geneSetID = geneset)

## Analysis of just DE genes
DE_df <- df %>% filter(FDR < 0.05) %>% filter(abs(logFC)>1)
DE_up <- DE_df %>% filter(logFC>0)
DE_down <- DE_df %>% filter(logFC<0)

## gProfiler
gprof_organism <- "mmusculus"
gostres <- gost(query = DE_up$ensembl_ID,
                organism = gprof_organism,
                user_threshold = 0.05,
                correction_method = "g_SCS")

## Interactive gProfiler plot
gostplot(gostres, capped = TRUE, interactive = TRUE)


## Use the next line to pick which terms you will highlight from the significant list
highlight_terms <- gostres$result$term_id   ## Highlight all
#highlight_terms <- gostres$result %>% slice_min(p_value,n = 25) %>% pull(term_id)   ## Highlight minimum p-value terms
highlight_terms <- gostres$result %>% filter(term_size < 1000) %>% slice_min(p_value,n = 25) %>% pull(term_id)   ## Highlight minimum p-value terms with term size < 1000 (gets rid of VERY general terms)

## Plot gProfiler for print
publish_gostplot(gostplot(gostres, capped = F, interactive = F), highlight_terms = highlight_terms, width = NA, height = NA, filename = NULL)


## If there are a lot of significant terms in gProfiler output, we can simplify with revigo
library(rrvgo)

## Calculate similarity between terms
simMatrix <- calculateSimMatrix(gostres$result$term_id,
                                orgdb="org.Mm.eg.db",
                                ont="BP",
                                method="Rel")

## Reduce terms and set names for plotting
scores <- setNames(-log10(gostres$result$p_value), gostres$result$term_id)
reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Mm.eg.db")

## Tree plot (grouped boxes)
treemapPlot(reducedTerms)


