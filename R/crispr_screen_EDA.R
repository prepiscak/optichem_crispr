# CRISPR analysis

# loading libraries ----
library("readr") # reading in files

library("dplyr") # data.frame manipulation
library("tidyr")
library("tibble") 
library("skimr") # summary for data.frame
# plots
library("ggplot2")
library("ggpubr")
library("pheatmap")

# Exploratory data analysis ----
DATASETS_DIR="./datasets/"
RESULTS_DIR="./results/"
dir.create(RESULTS_DIR)

# loading files ----
# TO-DO
# [ ] setup conda and renv!!!

# raw counts from plasmid, pretreatment, mock and tax
# [ ] sgRNA, gene and non-target coverage
# [ ] gini index etc.
# [ ] libA investigate genes with duplicated sequences!!! (Tceal1? or other genes from list?)

# all samples

raw_counts_all_file <- paste0(DATASETS_DIR, "total_all_tax_mock_counts/total_all_tax_mock.count.txt")
raw_counts_all <- readr::read_tsv(raw_counts_all_file)

count_summary_all_file <- paste0(DATASETS_DIR, "total_all_tax_mock_counts/total_all_tax_mock.countsummary.txt")
count_summary_all <- as.data.frame(readr::read_tsv(count_summary_all_file))

norm_counts_all_file <- paste0(DATASETS_DIR, "total_all_tax_mock_counts/total_all_tax_mock.count_normalized.txt")
norm_counts_all <- readr::read_tsv(norm_counts_all_file)

# mock and tax

# normalized (total) counts from from plasmid, pretreatment, mock and tax
# [ ] PCA, correlation plot
# [ ] filter and test for essential genes! (pathway enrichment): log2FC; removing NA padj (keep also < 2 sgRNAs?); perform only over-representation analysis; consider all detected genes as universe
raw_counts_file <- paste0(DATASETS_DIR, "total_tax_mock_counts/total_tax_mock.count.txt")
raw_counts <- readr::read_tsv(raw_counts_file)

# count_summary_file <- paste0(DATASETS_DIR, "total_tax_mock_counts/total_tax_mock.countsummary.txt")
# count_summary <- as.data.frame(readr::read_tsv(count_summary_file))

norm_counts_file <- paste0(DATASETS_DIR, "total_tax_mock_counts/total_tax_mock.count_normalized.txt")
norm_counts <- readr::read_tsv(norm_counts_file)

# re-normalized (total) counts only for mock and tax
# [ ] filter and run mageck test

# 


# information about sgRNAs and genes covered ----
# sgRNAs, genes
# overlap of not detected genes between condition, missing genes function (biological space)
# ignore NonTargetingControlGuideForMouse; especially for genes
# coverage of protein-coding and miRNAs (66405 sgRNAs) and non-target (1000) separately
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

# Final table for raw counts 

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
                total_sgRNAs=TotalsgRNAs) 

# descriptive analysis
# lower mapping 
# enough coverage for sgRNAs (not underrepresented?)
# tick some boxes for QC as on the page https://sourceforge.net/p/mageck/wiki/Home/#056
# summary statistics in conditions only: mock and tax
count_summary_all_clean_cond_only <- count_summary_all_clean %>%
  dplyr::filter(condition != "plasmid" & condition != "pretrt")

skimr::skim(count_summary_all_clean_cond_only)

readr::write_tsv(count_summary_all_clean, 
          path = paste0(RESULTS_DIR, "count_summary_all.tsv"))

# normalized counts ----
norm_counts_all_log2 <- norm_counts_all %>%
  dplyr::select(-Gene) %>%
  dplyr::mutate_at(vars(-matches("sgRNA")), function(x) log2(x+1)) %>%
  tibble::column_to_rownames(., var="sgRNA") 

labels_all_levels <- colnames(norm_counts_all_log2)
norm_counts_all_log2_gather <- norm_counts_all_log2 %>%
  tidyr::gather(., key="sample_label", value="log2counts") %>%
  dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1", sample_label)) %>%
  dplyr::mutate(condition = factor(condition, levels=c("plasmid", "pretrt", "mock", "tax"))) %>%
  dplyr::mutate(sample_label = factor(sample_label, levels = labels_all_levels))


norm_counts_all_dens_plot <- ggpubr::ggdensity(norm_counts_all_log2_gather, x = "log2counts",
                                           #add = "median", 
                                           #add = "mean",
                                           rug = FALSE,
                                           alpha = 0.1,
                                           fill = "sample_label",
                                           facet.by = "condition",
                                           panel.labs = list(condition = c("plasmid", "pre-treatment", "mock", "taxane")),
                                           xlab = "log2(counts+1)") +
  theme(legend.position="none")

# norm_counts_dens_plot <- ggplot(norm_counts_all_log2_gather, aes(x=log2counts, fill=sample_label)) +
#   geom_density(alpha=.2) + 
#   facet_grid(condition ~ .) + 
#   theme_bw() + 
#   theme(#legend.position="none",
#     #legend.text = element_text(size=10, face="bold"),
#     #text = element_text(size=8),
#     plot.title = element_text(hjust = 0.5),
#     #axis.title.y =  element_blank(),
#     #axis.title.x = element_blank(),
#     #axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
#     axis.text.y = element_text(face = "bold")) #axis.text.y = element_text(face = "bold", size = 6)

norm_counts_all_boxplot <- ggpubr::ggboxplot(norm_counts_all_log2_gather, x = "sample_label", y = "log2counts",
               fill = "condition") +
  ylab("log2(counts+1)") + 
  theme(legend.text = element_text(size=12),
    axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x =  element_blank())
  
ggsave(norm_counts_all_dens_plot, 
       filename = paste0(RESULTS_DIR, "norm_counts_all_dens_plot.pdf"))

ggsave(norm_counts_all_boxplot, 
       filename = paste0(RESULTS_DIR, "norm_counts_all_boxplot.pdf"))

save(count_summary_all_clean,
     norm_counts_all_log2_gather,
     norm_counts_all_dens_plot, 
     norm_counts_all_boxplot, 
     file = paste0(RESULTS_DIR, "crispr_screen_EDA.RData"))

# norm_counts_boxplot <- ggplot(norm_counts_all_log2_gather, aes(x=sample_label, y=log2counts, fill=condition)) +
#   geom_boxplot() + theme_bw() #+ #facet_grid(condition ~ .) + 
# theme(legend.position="none") # Remove legend

# non-target

# norm_counts_nontarget_log2 <- norm_counts %>%
#   dplyr::filter(grepl(pattern = "NonTargetingControlGuideForMouse", Gene)) %>%
#   dplyr::select(-Gene) %>%
#   tibble::column_to_rownames(., var="sgRNA") %>%
#   mutate_all(., function(x) log2(x+1))
# 
# norm_counts_nontarget_log2_gather <- norm_counts_nontarget_log2 %>%
#   tidyr::gather(., key="sample_label", value="log2counts") %>%
#   dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1", sample_label)) %>%
#   dplyr::mutate(condition = factor(condition, levels=c("plasmid", "pretrt", "mock", "tax")))
# 
# norm_counts_nontarget_dens_plot <- ggplot(norm_counts_nontarget_log2_gather, aes(x=log2counts, fill=sample_label)) +
#   geom_density(alpha=.2) + facet_grid(condition ~ .) + theme_bw() + 
#   theme(legend.position="none") # Remove legend
# 
# norm_counts_nontarget_boxplot <- ggplot(norm_counts_nontarget_log2_gather, aes(x=sample_label, y=log2counts, fill=condition)) +
#   geom_boxplot() + theme_bw() #+ #facet_grid(condition ~ .) + 
# theme(legend.position="none") # Remove legend


#################### Principal component analysis ################

ntop_genes=nrow(norm_counts_all_log2) # default 500
#rld_counts_ntop <- rnaSelectTopVarGenes_var(rld_counts, ntop = ntop_genes)
norm_counts_all_log2_ntop <- rnaSelectTopVarGenes(norm_counts_all_log2, ntop = 65959)

pca_all <- stats::prcomp(t(norm_counts_all_log2), center = TRUE)



# if less than 5 components needed plot otherwise skip; maybe update up to 5 components?
pca_all_data <- data.frame(PC1 = pca_all$x[, 1], PC2 = pca_all$x[, 2], sample_label=names(pca_all$x[, 1])) %>%
  dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1", sample_label)) %>%
  dplyr::mutate(batch)




pcaData <- plotPCA(rld_filt, intgroup=c("condition", "date_preparation"), returnData=TRUE, ntop=16821)
percentVar <- round(100 * attr(pcaData, "percentVar"))
PCA_individual_effect <- ggplot(pcaData, aes(PC1, PC2, color=condition, shape=date_preparation)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) #+ coord_fixed()

# Principal component analysis and correlation plots ----

pca <- prcomp(t(norm_counts_all_log2),center=TRUE)

percentVar <- pca$sdev^2/sum(pca$sdev^2)
#condition = run_details$Condition[run_details$Sample.names == colnames(slmat)]

d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample_label=names(pca$x[, 1])) %>%
  dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1", sample_label)) 
  
pca_all = ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "condition")) + 
  geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed() + 
  theme_bw()

#################### Correlation matrix ################
### Compute pairwise correlation values
correlation_breaks = seq(0, 1, 0.02)
legend_breaks = seq(0, 1, 0.2)
legend_breaks=NA
#correlation_breaks=NA
color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(correlation_breaks))

# Pearson correlation ----
norm_counts_all_log2_cor <- cor(norm_counts_all_log2, method = "pearson")    ## cor() is a base R function

S2S_pearson_plot <- pheatmap(norm_counts_all_log2_cor, color = color, border_color=NA, 
                             fontsize = 10, fontsize_row = 10, height=20, 
                             breaks = correlation_breaks,
                             legend_breaks = legend_breaks,
                             silent = TRUE,
                             main = "Sample-to-sample pearson correlation") # silent - do not draw the plot (useful when using the gtable output)



# Spearman correlation ----
norm_counts_all_log2__spearman_cor <- cor(norm_counts_all_log2, method = "spearman")    ## cor() is a base R function

S2S_spearman_plot <- pheatmap(norm_counts_all_log2__spearman_cor, color = color, border_color=NA, 
                              fontsize = 10, fontsize_row = 10, height=20, 
                              breaks = correlation_breaks,
                              legend_breaks = legend_breaks,
                              silent = TRUE,
                              main = "Sample-to-sample spearman correlation") # silent - do not draw the plot (useful when using the gtable output)

# saving data ----
ggsave(filename = paste0(RESULTS_DIR, "S2S_pearson_plot.pdf"), plot=S2S_pearson_plot)
ggsave(filename = paste0(RESULTS_DIR, "S2S_spearman_plot.pdf"), plot=S2S_spearman_plot)

# Saving S2S plots to RData object
save(count_summary_all_clean,
     norm_counts_all_log2,
     norm_counts_all_log2_gather,
     norm_counts_all_dens_plot, 
     norm_counts_all_boxplot, 
     S2S_pearson_plot,
     S2S_spearman_plot,
     file = paste0(RESULTS_DIR, "crispr_screen_EDA.RData"))


  
 















# test mageck-flute
library("MAGeCKFlute")

MapRatesView(count_summary)
IdentBarView(count_summary, x = "Label", y = "GiniIndex", 
             ylab = "Gini index", main = "Evenness of sgRNA reads")
count_summary$Missed = log10(count_summary$Zerocounts)
IdentBarView(count_summary, x = "Label", y = "Missed", fill = "#394E80",
             ylab = "Log10 missed gRNAs", main = "Missed sgRNAs")



libA_file <- "/media/prepisca/DataAnalysis1/optichem_crispr/reanalysis/DatasetsRaw/mock_tax/fastq/trimmed_fastq/Mouse_GeCKOv2_libA_09Mar2015_ForMageck.csv"
libA <- read_csv(libA_file, col_names = c("sgRNA", "sequence", "Gene"))

raw_counts_file <- "/media/prepisca/DataAnalysis1/optichem_crispr/reanalysis/DatasetsRaw/mock_tax/fastq/trimmed_fastq/total_tax_mock_counts/total_tax_mock_all_counts/total_tax_mock_all.count.txt"
raw_counts <- read_tsv(raw_counts_file)

norm_counts_file <- "/media/prepisca/DataAnalysis1/optichem_crispr/reanalysis/DatasetsRaw/mock_tax/fastq/trimmed_fastq/total_tax_mock_counts/total_tax_mock_all_counts/total_tax_mock_all.count_normalized.txt"
norm_counts <- read_tsv(norm_counts_file)

#norm_counts_file <- "/media/prepisca/DataAnalysis1/optichem_crispr/reanalysis/DatasetsRaw/mock_tax/fastq/trimmed_fastq/total_tax_mock_counts/total_tax_mock.count_normalized.txt"
#norm_counts <- read_tsv(norm_counts_file)

# Re-run normalization with mageck v0.5.6

# EDA

# library A exploration ----
# targets with identical sequences
# libB 62804 sgRNAs
# libA 67405; 22786 = 20611 (coding) + 1175 (miRNA) + 1000 (control non-target)
# 3sgRNAs per gene and 4 sgRNAs per miRNA:
# 61833 sgRNAs (coding) + 4700 (miRNA) + 1000 (control non-target)

libA_nontarget <- libA %>%
  dplyr::filter(grepl(pattern = "NonTargetingControlGuideForMouse", Gene))

libA_genes_mirna <- libA %>%
  dplyr::filter(!(Gene %in% libA_nontarget$Gene))

libA_genes <- libA_genes_mirna %>%
  dplyr::filter(!grepl(pattern = "mmu-", Gene))



# add t-test; wilcoxon test or similar for Tceal1 and also plot Tceal1 profiles
tceal1_exp <- norm_counts %>%
  dplyr::filter(Gene=="Tceal1") %>%
  dplyr::select(-Gene) %>%
  #tibble::column_to_rownames(., var="sgRNA") %>%
  dplyr::mutate_at(vars(-matches("sgRNA")), function(x) log2(x+1)) 
#%>% #tibble::rownames_to_column(., var = "sgRNA")

tceal1_exp_gather <- tceal1_exp %>%
  tidyr::gather(., key="sample_label", value="log2counts", -sgRNA) %>%
  dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1", sample_label)) %>%
  dplyr::filter(condition != "plasmid" & condition != "pretrt") %>%
  dplyr::mutate(condition = factor(condition, levels=c("mock", "tax")))

tceal1_sgrnas <- ggplot(tceal1_exp_gather, aes(x=sgRNA, y=log2counts, fill=condition)) +
  geom_boxplot() + 
  geom_line(aes(group = interaction(index, variable)),
            alpha = 0.5, colour = "darkgrey")
  theme_bw() 
#+ facet_grid(sgRNA ~ ., scales = "free_y")  #+ 
theme(legend.position="none") # Remove legend

tceal1_sgrnas <- ggplot(tceal1_exp_gather, aes(x=sgRNA, y=log2counts, fill=condition)) +
  geom_boxplot() + 
theme_bw() 

tceal1_sgrnas <- ggpubr::ggboxplot(tceal1_exp_gather, x = "sgRNA", y = "log2counts",
                                    fill = "condition") +
  ylab("log2(counts+1)") + 
  theme(legend.text = element_text(size=12)) #axis.title.x =  element_blank()

# ggpubr::ggpaired(tceal1_exp_gather, x = "condition", y = "log2counts",
#                  color = "condition", line.color = "gray", line.size = 0.4,
#                  palette = "npg", facet.by="sgRNA")

  
# Extracting only mock and tax ----  
# filtering lowly expressed sgRNAs
# keep >= 100 norm. counts and > 1 sgRNA
# 1. keep nontarget and do above filtering and use --control-sgrna for generating the null distribution of RRA
# 2. remove nontarget and do not use --control-sgrna    
norm_counts_gather <- norm_counts %>%
    tidyr::gather(., key="sample_label", value="norm_counts", -sgRNA, -Gene) %>%
    dplyr::mutate(condition=gsub(pattern="(.+)(_VH.+)", replacement = "\\1", sample_label)) # generating conditions

# 1. keep nontarget and do above filtering and use --control-sgrna for generating the null distribution of RRA     
norm_counts_grouped <- norm_counts_gather %>%
 dplyr::group_by(condition, sgRNA) %>%
 dplyr::summarise(gene=unique(Gene),
                  mean=mean(norm_counts)) %>% # using mean expression (potential outliers)
 tidyr::spread(., condition, mean)

# norm_counts_grouped <- norm_counts_gather %>%
#   dplyr::group_by(condition, sgRNA) %>%
#   dplyr::summarise(gene=unique(Gene),
#                    median=median(norm_counts)) %>% # using median instead of mean
#   tidyr::spread(., condition, median)

min_expression=100
min_sgrnas=2
cat("sgRNAs before expression filtering:", nrow(norm_counts_grouped), "\n")

# filtering based on expression ----
norm_counts_grouped_expfilt <- norm_counts_grouped %>%
  dplyr::filter( (mock >= min_expression) | (tax >= min_expression))
    
cat("sgRNAs after expression filtering:", nrow(norm_counts_grouped_expfilt), "\n")

# filtering based on number of sgrnas ----
norm_counts_grouped_filt <- norm_counts_grouped_expfilt %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(sgrnas=paste(unique(sgRNA), collapse = ";"),
                   number_sgrnas=n()) %>%
  #dplyr::filter(number_sgrnas >= min_sgrnas) # removing control non-target
  dplyr::filter(number_sgrnas >= min_sgrnas | grepl(pattern = "NonTargetingControlGuideForMouse", gene)) # keeping control non-target

keep_sgrnas <- norm_counts_grouped_filt %>%
  tidyr::separate_rows(., sgrnas, sep = ";") %>% # separate delimited values into respective rows
  dplyr::pull(sgrnas)

cat("sgRNAs after expression and sgRNA filtering:", length(keep_sgrnas), "\n")
cat("genes after expression and sgRNA filtering:", nrow(norm_counts_grouped_filt), "\n")
cat("  from which control nontarget:", sum(grepl(pattern = "NonTargetingControlGuideForMouse", norm_counts_grouped_filt$gene)), "\n")

# 1. keep nontarget and do above filtering and use --control-sgrna for generating the null distribution of RRA  
# filtering only on expression
norm_counts_filt_expression <- norm_counts %>%
  dplyr::filter(sgRNA %in% norm_counts_grouped_expfilt$sgRNA)

# filtering on expression and number of sgRNAs and keeping control non-target
norm_counts_filt_exp_sgrnas <- norm_counts %>%
  dplyr::filter(sgRNA %in% keep_sgrnas)

# filtering on expression and number of sgRNAs and removing control non-target
norm_counts_filt_exp_sgrnas_rm_nontarget <- norm_counts %>%
  dplyr::filter(sgRNA %in% keep_sgrnas) %>%
  dplyr::filter(!grepl("NonTargetingControlGuideForMouse", Gene))

# filtering raw counts and re-normalazing in mageck test!
#  on expression and number of sgRNAs and removing control non-target
raw_counts_filt_exp_sgrnas_rm_nontarget <- raw_counts %>%
  dplyr::filter(sgRNA %in% keep_sgrnas) %>%
  dplyr::filter(!grepl("NonTargetingControlGuideForMouse", Gene))

# filtered only based on expression
raw_counts_filt_exp_rm_nontarget <- raw_counts %>%
  dplyr::filter(sgRNA %in% norm_counts_grouped_expfilt$sgRNA) %>%
  dplyr::filter(!grepl("NonTargetingControlGuideForMouse", Gene))

readr::write_tsv(norm_counts_filt_expression, 
                 path = paste0(RESULTS_DIR, "norm_counts_filt_expression.tsv"))

readr::write_tsv(norm_counts_filt_exp_sgrnas, 
                 path = paste0(RESULTS_DIR, "norm_counts_filt_exp_sgrnas.tsv"))

readr::write_tsv(norm_counts_filt_exp_sgrnas_rm_nontarget, 
                 path = paste0(RESULTS_DIR, "norm_counts_filt_exp_sgrnas_rm_nontarget.tsv"))

readr::write_tsv(raw_counts_filt_exp_sgrnas_rm_nontarget, 
                 path = paste0(RESULTS_DIR, "raw_counts_filt_exp_sgrnas_rm_nontarget.tsv"))

readr::write_tsv(raw_counts_filt_exp_rm_nontarget, 
                 path = paste0(RESULTS_DIR, "raw_counts_filt_exp_rm_nontarget.tsv"))



# plot distribution of non-targets after filtering ----
# same as below, but for genes rather than control non-target

newly_norm_counts <- readr::read_tsv("/media/prepisca/DataAnalysis1/optichem_crispr/optichem_crispr_tceal1/results/tax_mock_renorm/tax_mock_renorm.normalized.txt")

# plot distribution of non-targets after filtering ----
norm_counts_nontarget_exp_filt_log2 <- norm_counts_filt_expression %>%
  dplyr::filter(!grepl(pattern = "NonTargetingControlGuideForMouse", Gene)) %>%
  dplyr::select(-Gene) %>%
  tibble::column_to_rownames(., var="sgRNA") %>%
  mutate_all(., function(x) log2(x+1))


pca <- prcomp(t(norm_counts_nontarget_exp_filt_log2),center=TRUE)

percentVar <- pca$sdev^2/sum(pca$sdev^2)
#condition = run_details$Condition[run_details$Sample.names == colnames(slmat)]

d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample_label=names(pca$x[, 1])) %>%
  dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1", sample_label)) 

pca_all = ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "condition")) + 
  geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
  ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + coord_fixed()


norm_counts_nontarget_exp_filt_log2_gather <- norm_counts_nontarget_exp_filt_log2 %>%
  tidyr::gather(., key="sample_label", value="log2counts") %>%
  dplyr::mutate(condition=gsub(pattern = "(.+)(_VH.+)", replacement = "\\1", sample_label)) %>%
  dplyr::mutate(condition = factor(condition, levels=c("plasmid", "pretrt", "mock", "tax")))

norm_counts_nontarget_exp_filt_dens_plot <- ggplot(norm_counts_nontarget_exp_filt_log2_gather, aes(x=log2counts, fill=sample_label)) +
  geom_density(alpha=.2) + facet_grid(condition ~ .) + theme_bw() #+ 
  theme(legend.position="none") # Remove legend

norm_counts_nontarget_exp_filt_boxplot <- ggplot(norm_counts_nontarget_exp_filt_log2_gather, aes(x=sample_label, y=log2counts, fill=condition)) +
  geom_boxplot() + theme_bw() #+ #facet_grid(condition ~ .) + 
theme(legend.position="none") # Remove legend

