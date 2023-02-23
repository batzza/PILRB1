load("data/rnaseq/processed/Processed_RNAseq_w_CPM.Rdata")

library(dplyr)

## Load gene lists
ptc_genes <- read.csv("data/rnaseq/external/enriched_gene_groups/Reactome_PTC_mouse.txt", sep = '\t')
bipolar_genes <- read.csv("data/rnaseq/external/enriched_gene_groups/HPA/mouse_orthologs/Bipolar_enriched_mouse.txt", sep = '\t') %>% filter(Enriched_type == "Bipolar cells")
bipolar_only_genes <- bipolar_genes %>% filter(Enriched_type_count==1)
cone_genes <- read.csv("data/rnaseq/external/enriched_gene_groups/HPA/mouse_orthologs/Cone_enriched_mouse.txt", sep = '\t') %>% filter(Enriched_type == "Cone photoreceptor cells")
cone_only_genes <- cone_genes %>% filter(Enriched_type_count==1)
horizontal_genes <- read.csv("data/rnaseq/external/enriched_gene_groups/HPA/mouse_orthologs/Horizontal_enriched_mouse.txt", sep = '\t') %>% filter(Enriched_type == "Horizontal cells")
horizontal_only_genes <- cone_genes %>% filter(Enriched_type_count==1)
rod_genes <- read.csv("data/rnaseq/external/enriched_gene_groups/HPA/mouse_orthologs/Rod_enriched_mouse.txt", sep = '\t') %>% filter(Enriched_type == "Rod photoreceptor cells")
rod_only_genes <- rod_genes %>% filter(Enriched_type_count==1)

GPCR_genes <- read.csv("data/rnaseq/external/enriched_gene_groups/NHLBI_GPCR_mouse_list.txt", sep = '\t')
Ribbon_genes <- read.csv("data/rnaseq/external/enriched_gene_groups/Ribbon_synapse_mouse.txt", sep = '\t')



summarize_counts <- function(all_table, sig_table, gene_list, list_name, age){
  n_exp <- all_table %>% filter(ensembl_ID %in% gene_list$Ensembl) %>% nrow()
  n_not_expressed <- nrow(gene_list) - n_exp
  n_DE_up_in_KO <- sig_table %>% filter(ensembl_ID %in% gene_list$Ensembl) %>% filter(logFC > 0) %>% nrow()
  n_DE_down_in_KO <- sig_table %>% filter(ensembl_ID %in% gene_list$Ensembl) %>% filter(logFC < 0) %>% nrow()
  n_no_change <- n_exp - n_DE_up_in_KO - n_DE_down_in_KO
  
  out <- data.frame(age = age,
             gene_list  = list_name,
             n_genes = nrow(gene_list),
             n_exp = n_exp,
             n_not_exp = n_not_expressed,
             n_DE_up_in_KO = n_DE_up_in_KO,
             n_DE_down_in_KO = n_DE_down_in_KO,
             n_no_change = n_no_change
  )
  
  return(out)
}

#P28
count_by_list <- summarize_counts(all_table = DE$P28_WTvKO, sig_table = DE$P28_WTvKO_sig, gene_list = ptc_genes, list_name = "PTC", age = "P28")
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$P28_WTvKO, sig_table = DE$P28_WTvKO_sig, gene_list = bipolar_genes, list_name = "Bipolar_enriched", age = "P28"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$P28_WTvKO, sig_table = DE$P28_WTvKO_sig, gene_list = cone_genes, list_name = "Cone_enriched", age = "P28"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$P28_WTvKO, sig_table = DE$P28_WTvKO_sig, gene_list = horizontal_genes, list_name = "Horizontal_enriched", age = "P28"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$P28_WTvKO, sig_table = DE$P28_WTvKO_sig, gene_list = rod_genes, list_name = "Rod_enriched", age = "P28"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$P28_WTvKO, sig_table = DE$P28_WTvKO_sig, gene_list = GPCR_genes, list_name = "GPCR", age = "P28"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$P28_WTvKO, sig_table = DE$P28_WTvKO_sig, gene_list = Ribbon_genes, list_name = "Ribbon", age = "P28"))

#M12
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$M12_WTvKO, sig_table = DE$M12_WTvKO_sig, gene_list = ptc_genes, list_name = "PTC", age = "M12"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$M12_WTvKO, sig_table = DE$M12_WTvKO_sig, gene_list = bipolar_genes, list_name = "Bipolar_enriched", age = "M12"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$M12_WTvKO, sig_table = DE$M12_WTvKO_sig, gene_list = cone_genes, list_name = "Cone_enriched", age = "M12"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$M12_WTvKO, sig_table = DE$M12_WTvKO_sig, gene_list = horizontal_genes, list_name = "Horizontal_enriched", age = "M12"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$M12_WTvKO, sig_table = DE$M12_WTvKO_sig, gene_list = rod_genes, list_name = "Rod_enriched", age = "M12"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$M12_WTvKO, sig_table = DE$M12_WTvKO_sig, gene_list = GPCR_genes, list_name = "GPCR", age = "M12"))
count_by_list <- rbind(count_by_list, summarize_counts(all_table = DE$M12_WTvKO, sig_table = DE$M12_WTvKO_sig, gene_list = Ribbon_genes, list_name = "Ribbon", age = "M12"))

count_by_list


ptc_genes %>% filter(Ensembl %in% DE$M12_WTvKO_sig$ensembl_ID)
DE$P28_WTvKO %>% filter(ensembl_ID %in% bipolar_genes$Ensembl) %>% filter(FDR<0.05) %>% filter(abs(logFC) > 0)
DE$M12_WTvKO %>% filter(ensembl_ID %in% Ribbon_genes$Ensembl) %>% filter(FDR<0.05) %>% filter(abs(logFC) > 1)


DE$P28_WTvKO %>% filter(ensembl_ID %in% c("ENSMUSG00000040554"))
DE$M12_WTvKO %>% filter(ensembl_ID %in% c("ENSMUSG00000040554"))
