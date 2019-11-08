# CRISPR analysis

# A. crispr_screen_EDA.R
# B. crispr_screen_filtering.R
# C. crispr_screen_results.R
# ADD final results table of top N genes (FDR < 0.25)
# Methods
# Results description
# GSEA based on rank? (with without FDR filtering)

#################### Loading libraries ################
library("readr") # reading in files

library("dplyr") # data.frame manipulation
library("tidyr")
library("tibble") 

# plots
library("ggplot2")
library("ggrepel")
library("ggpubr")
library("pheatmap")
library("kableExtra") # to save data.frame as pdf

#################### Loading datasets ################
DATASETS_DIR="./datasets/"
RESULTS_DIR="./results/crispr_screen_results/"
dir.create(RESULTS_DIR, recursive=TRUE)

# gene summary file
rra_gene_summary_file <- paste0(DATASETS_DIR, "tax_mock_test/tax_mock_test.gene_summary.txt")
rra_gene_summary <- readr::read_tsv(rra_gene_summary_file, col_types="cidddiiddddiid")

# gene summary file
rra_sgrna_summary_file <- paste0(DATASETS_DIR, "tax_mock_test/tax_mock_test.sgrna_summary.txt")
rra_sgrna_summary <- readr::read_tsv(rra_sgrna_summary_file, col_types="ccccddddddddddc")

# normalized counts - for checking
rra_norm_counts_file <- paste0(DATASETS_DIR, "tax_mock_test/tax_mock_test.normalized.txt")
rra_norm_counts <- readr::read_tsv(rra_norm_counts_file)

#################### Exploratory data analysis ################

# extracting significant results in negative selection
padj_cutoff <- 0.25
rra_gene_summary_signif <- rra_gene_summary %>%
  dplyr::select(id, num, `neg|goodsgrna`, `neg|lfc`, `neg|fdr`) %>% #`neg|p-value`, 
  dplyr::rename(`gene symbol` = id,
                `detected sgRNAs` = num,
                `good sgRNAs` = `neg|goodsgrna`,
                `log2FC` = `neg|lfc`,
                `padj` = `neg|fdr`) %>%
  dplyr::filter(!is.na(padj) & padj < padj_cutoff)

# add extra information about gene? gene description?
kableExtra::kable(x = rra_gene_summary_signif,
                  format = "latex", 
                  digits=c(0,0,0,1,4),
                  booktabs = TRUE, 
                  linesep="") %>% #add linespace after each condition e.g. linesep=c("", "\\addlinespace"); or completely remove linesep=""
  column_spec(c(2,3), width = "4em") %>%
  save_kable(paste0(RESULTS_DIR, "rra_gene_summary_signif.pdf"))


#################### sgRNA counts ################

# log2 transformation
rra_signif_norm_counts_log2 <- rra_norm_counts %>%
  dplyr::filter(Gene %in% rra_gene_summary_signif$`gene symbol`) %>%
  dplyr::mutate_at(vars(-matches("sgRNA|Gene")), function(x) log2(x+1)) 
  
# keeping original sgRNA names
# rra_signif_norm_counts_log2_gather <- rra_signif_norm_counts_log2 %>%
#   tidyr::gather(., key="sample_label", value="log2counts", -sgRNA, -Gene) %>%
#   dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1", sample_label)) %>%
#   dplyr::mutate(condition = factor(condition, levels=c("mock", "tax")))

# renaming sgRNAs
# easier way of renaming?!
rra_signif_norm_counts_log2_gather <- rra_signif_norm_counts_log2 %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(sgRNA = paste("sgRNA", seq(n()), sep = "", collapse = ";"),
                   mock_VH1_IL = paste(mock_VH1_IL, collapse = ";"),
                   mock_VH3_I2L2R = paste(mock_VH3_I2L2R, collapse = ";"),
                   mock_VH4_II2R = paste(mock_VH4_II2R, collapse = ";"),
                   mock_VH5_IIL = paste(mock_VH5_IIL, collapse = ";"),
                   mock_VH6_IIR = paste(mock_VH6_IIR, collapse = ";"),
                   mock_VH7_IIA2R = paste(mock_VH7_IIA2R, collapse = ";"),
                   mock_VH8_IIALR = paste(mock_VH8_IIALR, collapse = ";"),
                   mock_VH9_IIAL = paste(mock_VH9_IIAL, collapse = ";"),
                   mock_VH10_I2L = paste(mock_VH10_I2L, collapse = ";"),
                   tax_VH15_IBR = paste(tax_VH15_IBR, collapse = ";"),
                   tax_VH16_IBL = paste(tax_VH16_IBL, collapse = ";"),
                   tax_VH17_IB2L = paste(tax_VH17_IB2L, collapse = ";"),
                   tax_VH23_IBLR = paste(tax_VH23_IBLR, collapse = ";"),
                   tax_VH24_IB2R = paste(tax_VH24_IB2R, collapse = ";")) %>%
  tidyr::separate_rows(., -Gene, sep=";") %>%
  tidyr::gather(., key="sample_label", value="log2counts", -sgRNA, -Gene) %>%
  dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1", sample_label),
                log2counts = as.numeric(log2counts)) %>%
  dplyr::mutate(condition = factor(condition, levels=c("mock", "tax"))) %>%
  dplyr::mutate(Gene = factor(Gene, levels = rra_gene_summary_signif$`gene symbol`))

# plot heatmap of significant
rra_signif_norm_counts_log2_heatmap <- rra_signif_norm_counts_log2 %>%
  dplyr::select(-Gene) %>%
  tibble::column_to_rownames(var="sgRNA") %>%
  pheatmap::pheatmap(.,
                     scale="none",
                     cluster_rows = FALSE,
                     cluster_cols = FALSE)


# sgRNA counts for significant genes
# potentially add t-test; wilcox test
rra_all_signif_norm_counts_log2_boxplot <- ggpubr::ggboxplot(rra_signif_norm_counts_log2_gather, 
                                             x = "sgRNA", 
                                             y = "log2counts",
                                             fill = "condition",
                                             facet.by = "Gene") +
  ylab("log2(counts+1)") + 
  theme(legend.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x =  element_blank())

ggsave(rra_all_signif_norm_counts_log2_boxplot, 
       width = 12,
       height = 11.5,
       #units = "in",
       filename = paste0(RESULTS_DIR, "rra_gene_boxplot_all_signif.pdf"))

selected_signif <- c("Tceal1", "Cul9", "Defa25", "Mettl10", "Ccl21a", "Wdr72", "Hist1h2bc")
rra_signif_norm_counts_log2_boxplot <- ggpubr::ggboxplot(dplyr::filter(rra_signif_norm_counts_log2_gather, Gene %in% selected_signif), 
                                                             x = "sgRNA", 
                                                             y = "log2counts",
                                                             fill = "condition",
                                                             facet.by = "Gene") +
  ylab("log2(counts+1)") + 
  theme(legend.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x =  element_blank())

ggsave(rra_signif_norm_counts_log2_boxplot, filename = paste0(RESULTS_DIR, "rra_gene_boxplot_signif.pdf"))

# Calculate sgRNA abundance
sgrna_abundance <- rra_norm_counts %>%
  dplyr::select(-Gene) %>%
  tibble::column_to_rownames(., var="sgRNA") 

total_sgrnas_condition <- apply(sgrna_abundance, MARGIN = 2, sum)
sgrna_abundance_freq <- sweep(sgrna_abundance, MARGIN = 2, STATS = total_sgrnas_condition, FUN = "/")
sgrna_abundance_freq_perc <- sgrna_abundance_freq*100

sgrna2gene_mapping <- rra_norm_counts %>%
  dplyr::select(sgRNA, Gene)
sgrna_abundance_freq_perc_gene <- sgrna_abundance_freq_perc %>%
  tibble::rownames_to_column(., var="sgRNA") %>%
  left_join(., sgrna2gene_mapping, by = "sgRNA") %>%
  dplyr::select(Gene, sgRNA, everything())

save(rra_gene_summary,
     rra_gene_summary_signif,
     rra_all_signif_norm_counts_log2_boxplot,
     rra_signif_norm_counts_log2_boxplot,
     file = paste0(RESULTS_DIR, "crispr_screen_results.RData"))