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
library("readxl")

library("dplyr") # data.frame manipulation
library("tidyr")
library("tibble") 
library("skimr") # summary for data.frame

# plots
library("ggplot2")
library("ggrepel")
library("ggpubr")
library("pheatmap")
library("kableExtra") # to save data.frame as pdf

#################### Loading datasets ################
DATASETS_DIR="./datasets/"
RESULTS_DIR="./results/crispr_screen_EDA/"
dir.create(RESULTS_DIR, recursive=TRUE)

# all samples metadata
metadata_all_file <- paste0(DATASETS_DIR, "crispr_invivo_metadata.xlsx")
metadata_all <- readxl::read_xlsx(metadata_all_file)

# all samples (plasmid, pretreatment, mock and tax) raw counts
raw_counts_all_file <- paste0(DATASETS_DIR, "total_all_tax_mock_counts/total_all_tax_mock.count.txt")
raw_counts_all <- readr::read_tsv(raw_counts_all_file)

# all samples total normalized counts
norm_counts_all_file <- paste0(DATASETS_DIR, "total_all_tax_mock_counts/total_all_tax_mock.count_normalized.txt")
norm_counts_all <- readr::read_tsv(norm_counts_all_file)

# all samples count summary 
count_summary_all_file <- paste0(DATASETS_DIR, "total_all_tax_mock_counts/total_all_tax_mock.countsummary.txt")
count_summary_all <- as.data.frame(readr::read_tsv(count_summary_all_file))

# mock and tax
# raw counts from mock and tax
#raw_counts_tax_file <- paste0(DATASETS_DIR, "total_tax_mock_counts/total_tax_mock.count.txt")
#raw_counts_tax <- readr::read_tsv(raw_counts_tax_file)

# total normalized counts
norm_counts_tax_file <- paste0(DATASETS_DIR, "total_tax_mock_counts/total_tax_mock.count_normalized.txt")
norm_counts_tax <- readr::read_tsv(norm_counts_tax_file)

#################### Exploratory data analysis ################
# [x] examining sgRNA and genes count distribution
# [ ]

#################### sgRNA, genes coverage ################
# information about sgRNAs and genes covered
# coverage of protein-coding and miRNAs (66405 sgRNAs) and non-target (1000)

# raw counts ----
raw_counts_all_coverage <- raw_counts_all %>%
  tidyr::gather(., key="sample_label", value="counts", -sgRNA, -Gene) %>%
  dplyr::filter(counts != 0) # filtering zero-counts

total_sgrnas=length(unique(raw_counts_all$sgRNA)) #65959
raw_counts_all_coverage_sgRNAs <- raw_counts_all_coverage %>%
  dplyr::group_by(sample_label) %>%
  dplyr::summarise(detected_sgrnas=n()) %>%
  dplyr::mutate(detected_sgrnas_percentage=(detected_sgrnas/total_sgrnas)*100)

# break-down into protein-coding, miRNAs, nontargetting
# total genes (theoretical): 22786 = 20611 (coding) + 1175 (miRNA) + 1000 (control non-target)
raw_counts_all_coverage_only_genes <- raw_counts_all_coverage %>%
  dplyr::filter(!grepl("NonTargetingControlGuideForMouse", Gene)) 

total_genes <- length(unique(raw_counts_all_coverage_only_genes$Gene)) # 21485
raw_counts_all_coverage_genes <- raw_counts_all_coverage_only_genes %>%
  dplyr::group_by(sample_label) %>%
  dplyr::distinct(Gene) %>%
  dplyr::summarise(detected_genes=n()) %>%
  dplyr::mutate(detected_genes_percentage=(detected_genes/total_genes)*100) %>%
  dplyr::mutate(total_genes=total_genes)

# control non-target coverage
# raw_counts_coverage_nontarget <- raw_counts_coverage %>%
#   dplyr::filter(grepl(pattern = "NonTargetingControlGuideForMouse", Gene))
# 
# raw_counts_coverage_nontarget_sgRNAs <- raw_counts_coverage_nontarget %>%
#   dplyr::group_by(sample_label) %>%
#   dplyr::summarise(detected_sgrnas=n()) %>%
#   dplyr::mutate(detected_sgrnas_perc=(detected_sgrnas/1000)*100)

# final table for raw counts 
count_summary_all_upd <- count_summary_all %>%
  dplyr::left_join(., raw_counts_all_coverage_sgRNAs, by = c("Label" = "sample_label")) %>%
  dplyr::left_join(., raw_counts_all_coverage_genes, by = c("Label" = "sample_label")) %>%
  dplyr::mutate(Reads_millions=(Reads/10^6),
                Mapped_millions=(Mapped/10^6),
                Mapped_percentage=Percentage*100)
  
# keeping only key information
count_summary_all_clean <- count_summary_all_upd %>%
  dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1" ,Label)) %>%
  dplyr::select(Label, condition, Reads_millions, Mapped_millions, Mapped_percentage, GiniIndex,
                TotalsgRNAs, detected_sgrnas, detected_sgrnas_percentage, 
                total_genes, detected_genes, detected_genes_percentage) %>%
  dplyr::rename(sample_label=Label,
                reads_millions=Reads_millions,
                mapped_millions=Mapped_millions,
                mapped_percentage=Mapped_percentage,
                gini_index=GiniIndex,
                total_sgRNAs=TotalsgRNAs) %>%
  dplyr::rename(`sample label` = sample_label,
                `reads [millions]` = reads_millions,
                `mapped [millions]` = mapped_millions,
                `mapped [%]` = mapped_percentage,
                `gini index` = gini_index,
                `total sgRNAs` = total_sgRNAs,
                `detected sgRNAs` = detected_sgrnas,
                `detected sgRNAs [%]` = detected_sgrnas_percentage,
                `total genes` = total_genes,
                `detected genes` = detected_genes,
                `detected genes [%]` = detected_genes_percentage) # rename for exporting

# summary statistics in conditions only: mock and tax
count_summary_all_clean_cond_only <- count_summary_all_clean %>%
  dplyr::filter(condition != "plasmid" & condition != "pretrt")

skimr::skim(count_summary_all_clean)
skimr::skim(count_summary_all_clean_cond_only)

# saving table
kableExtra::kable(x = count_summary_all_clean,
                  format = "latex", 
                  digits=1,
                  booktabs = TRUE, 
                  linesep="") %>% #add linespace after each condition e.g. linesep=c("", "\\addlinespace"); or completely remove linesep=""
  kable_styling(latex_options = c("scale_down")) %>% #" striped", 
  column_spec(c(3,4,7,8), width = "4em") %>%
  column_spec(c(9), width = "5.3em") %>%
  column_spec(c(12), width = "5.2em") %>%
  column_spec(c(5,11), width = "4em") %>%
  column_spec(c(6,10), width = "3em") %>%
  save_kable(paste0(RESULTS_DIR, "count_summary_all.pdf"))

# format="html", grey lines
kableExtra::kable(x = count_summary_all_clean, 
                  format = "html", digits=1,
                  booktabs = T) %>% #add linespace after each condition e.g. linesep=c("", "\\addlinespace"); or completely remove linesep=""
  kable_styling(bootstrap_options = "condensed") %>% #" striped", 
  save_kable(paste0(RESULTS_DIR, "count_summary_all2.pdf"))

readr::write_tsv(count_summary_all_clean, 
          path = paste0(RESULTS_DIR, "count_summary_all.tsv"))

# normalized counts ----
# creating matrix for PCA and correlation plot
norm_counts_all_log2 <- norm_counts_all %>%
  dplyr::select(-Gene) %>%
  dplyr::mutate_at(vars(-matches("sgRNA")), function(x) log2(x+1)) %>%
  tibble::column_to_rownames(., var="sgRNA") 

# keeping gene names to highlight control non-target
norm_counts_all_log2_gene <- norm_counts_all %>%
  dplyr::mutate_at(vars(-matches("sgRNA|Gene")), function(x) log2(x+1))

labels_all_levels <- colnames(norm_counts_all_log2) # sample_labels
norm_counts_all_log2_gene_gather <- norm_counts_all_log2_gene %>%
  tidyr::gather(., key="sample_label", value="log2counts", -sgRNA, -Gene) %>%
  dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1", sample_label)) %>%
  dplyr::mutate(condition = factor(condition, levels=c("plasmid", "pretrt", "mock", "tax"))) %>%
  dplyr::mutate(sample_label = factor(sample_label, levels = labels_all_levels))

# add highlight non-control target distribution! also in the boxplot below!
norm_counts_all_dens_plot <- ggpubr::ggdensity(norm_counts_all_log2_gene_gather, 
                                               x = "log2counts",
                                               #add = "median", 
                                               #add = "mean",
                                               rug = FALSE,
                                               alpha = 0.1,
                                               fill = "sample_label",
                                               facet.by = "condition",
                                               panel.labs = list(condition = c("plasmid", "pre-treatment", "mock", "taxane")),
                                               xlab = "log2(counts+1)") +
  #geom_density(data=dplyr::filter(norm_counts_all_log2_gene_gather, grepl(pattern = "NonTargetingControlGuideForMouse", x=Gene)), 
  #             aes(x = log2counts), colour="red") + # adding non-target; but for mock and tax it is summarized across condition
  theme(legend.position="none")

norm_counts_all_boxplot <- ggpubr::ggboxplot(norm_counts_all_log2_gene_gather, 
                                             x = "sample_label", 
                                             y = "log2counts",
                                             fill = "condition") +
  ylab("log2(counts+1)") + 
  theme(legend.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x =  element_blank())

norm_counts_all_violin <- ggpubr::ggviolin(norm_counts_all_log2_gene_gather, 
                                           x = "sample_label", 
                                           y = "log2counts",
                                           fill = "condition") +
  ylab("log2(counts+1)") + 
  theme(legend.text = element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x =  element_blank())

ggsave(norm_counts_all_dens_plot, 
       filename = paste0(RESULTS_DIR, "norm_counts_all_dens_plot.pdf"))

ggsave(norm_counts_all_boxplot, 
       filename = paste0(RESULTS_DIR, "norm_counts_all_boxplot.pdf"))

ggsave(norm_counts_all_violin, 
       filename = paste0(RESULTS_DIR, "norm_counts_all_violin.pdf"))

#################### Principal component analysis ################

# all samples
pca_all <- stats::prcomp(t(norm_counts_all_log2),center=TRUE)

pca_all_percentVar <- pca_all$sdev^2/sum(pca_all$sdev^2)

pca_all_data <- data.frame(PC1 = pca_all$x[, 1], 
                           PC2 = pca_all$x[, 2], 
                           sample_label=names(pca_all$x[, 1])) %>%
  dplyr::left_join(., metadata_all, by="sample_label") %>%
  dplyr::mutate(condition = if_else(sample_label == "pretrt_VH2VH12VH22", "pretrt", condition)) %>% # populating condition for merged pre-treatment
  dplyr::mutate(condition = factor(condition, levels = c("plasmid", "pretrt", "mock", "tax")))

pca_all_plot <- ggscatter(data = pca_all_data, 
                          x = "PC1", y = "PC2", 
                          color = "condition", 
                          #label = "sample_label",
                          repel = TRUE,
                          size=5,
                          xlab = paste0("PC1: ", round(pca_all_percentVar[1] * 100), "% variance"),
                          ylab = paste0("PC2: ", round(pca_all_percentVar[2] * 100), "% variance")) 

# mock and tax 
norm_counts_tax_log2 <- norm_counts_tax %>%
  dplyr::select(-Gene) %>%
  dplyr::mutate_at(vars(-matches("sgRNA")), function(x) log2(x+1)) %>%
  tibble::column_to_rownames(., var="sgRNA") 

pca_tax <- stats::prcomp(t(norm_counts_tax_log2),center=TRUE)

pca_tax_percentVar <- pca_tax$sdev^2/sum(pca_tax$sdev^2)

pca_tax_data <- data.frame(PC1 = pca_tax$x[, 1], 
                           PC2 = pca_tax$x[, 2], 
                           sample_label=names(pca_tax$x[, 1])) %>%
  dplyr::left_join(., metadata_all, by="sample_label") 


pca_tax_plot <- ggscatter(data = pca_tax_data, 
                          x = "PC1", y = "PC2", 
                          color = "condition", 
                          #label = "sample_label",
                          repel = TRUE,
                          size=5,
                          xlab = paste0("PC1: ", round(pca_tax_percentVar[1] * 100), "% variance"),
                          ylab = paste0("PC2: ", round(pca_tax_percentVar[2] * 100), "% variance")) 

# saving plots
ggsave(filename = paste0(RESULTS_DIR, "pca_all_plot.pdf"), plot=pca_all_plot)
ggsave(filename = paste0(RESULTS_DIR, "pca_tax_plot.pdf"), plot=pca_tax_plot)

#################### Correlation matrix ################
### Compute pairwise correlation values
correlation_breaks = seq(0, 1, 0.02)
legend_breaks = seq(0, 1, 0.2)
legend_breaks=NA
#correlation_breaks=NA
color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(correlation_breaks))

# pearson correlation
norm_counts_all_log2_cor <- cor(norm_counts_all_log2, method = "pearson")    ## cor() is a base R function

S2S_pearson_plot <- pheatmap(norm_counts_all_log2_cor, color = color, border_color=NA, 
                             fontsize = 10, fontsize_row = 10, height=20, 
                             breaks = correlation_breaks,
                             legend_breaks = legend_breaks,
                             silent = TRUE,
                             main = "Sample-to-sample pearson correlation") # silent - do not draw the plot (useful when using the gtable output)

# spearman correlation
norm_counts_all_log2__spearman_cor <- cor(norm_counts_all_log2, method = "spearman")    ## cor() is a base R function

S2S_spearman_plot <- pheatmap(norm_counts_all_log2__spearman_cor, color = color, border_color=NA, 
                              fontsize = 10, fontsize_row = 10, height=20, 
                              breaks = correlation_breaks,
                              legend_breaks = legend_breaks,
                              silent = TRUE,
                              main = "Sample-to-sample spearman correlation") # silent - do not draw the plot (useful when using the gtable output)

# saving plots and data
ggsave(filename = paste0(RESULTS_DIR, "S2S_pearson_plot.pdf"), plot=S2S_pearson_plot)
ggsave(filename = paste0(RESULTS_DIR, "S2S_spearman_plot.pdf"), plot=S2S_spearman_plot)

save(count_summary_all_clean,
     norm_counts_all_log2,
     norm_counts_all_log2_gene_gather,
     norm_counts_all_dens_plot, 
     norm_counts_all_boxplot, 
     norm_counts_all_violin,
     norm_counts_tax_log2,
     pca_all_plot,
     pca_tax_plot,
     S2S_pearson_plot,
     S2S_spearman_plot,
     file = paste0(RESULTS_DIR, "crispr_screen_EDA.RData"))