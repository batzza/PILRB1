## 01-05-2023
## Downstream analysis of gene expression changes in PILRB1 KO mice
## And figure generation

## Load the upstream processed data
load("data/rnaseq/processed/Processed_RNAseq_w_CPM.Rdata")

## load libraries
library(gprofiler2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(ggplot2)

## Set focal comparison 
df <- DE$P28_WTvKO

## Get log2 fold change 
gene_list <- df$logFC

## name the vector
names(gene_list) <- df$ensembl_ID

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)


## Run clusterProfiler GSEA
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Mm.eg.db, 
             pAdjustMethod = "BH")


## Dotplot
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)

## Ridgeplot
ridgeplot(gse) + labs(x = "enrichment distribution")

## GSEA plot
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)


# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(names(gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = df[df$ensembl_ID %in% dedup_ids$ENSEMBL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$ENTREZID = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$logFC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$ENTREZID

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)


kegg_organism = "mmu"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
