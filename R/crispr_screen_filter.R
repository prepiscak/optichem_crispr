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

#################### Loading datasets ################
DATASETS_DIR="./datasets/"
RESULTS_DIR="./datasets/crispr_screen_filter/"
dir.create(RESULTS_DIR, recursive=TRUE)

# mock and tax
# raw counts from mock and tax
raw_counts_tax_file <- paste0(DATASETS_DIR, "total_tax_mock_counts/total_tax_mock.count.txt")
raw_counts_tax <- readr::read_tsv(raw_counts_tax_file)

# total normalized counts
norm_counts_tax_file <- paste0(DATASETS_DIR, "total_tax_mock_counts/total_tax_mock.count_normalized.txt")
norm_counts_tax <- readr::read_tsv(norm_counts_tax_file)

#################### sgRNA filtering ################
# filtering lowly expressed sgRNAs
# filtering genes with less than 2 sgRNAs detected
# keep if on average >= 100 norm. counts in a condition group (mock or tax) and after this keep only genes for which > 1 sgRNA
# filter from sgRNAs from raw_counts, remove control nontarget and re-normalize    
norm_counts_tax_gather <- norm_counts_tax %>%
    tidyr::gather(., key="sample_label", value="norm_counts", -sgRNA, -Gene) %>%
    dplyr::mutate(condition=gsub(pattern="(.+)(_VH.+)", replacement = "\\1", sample_label)) # generating conditions

# 1. keep nontarget and do above filtering and use --control-sgrna for generating the null distribution of RRA     
norm_counts_tax_grouped <- norm_counts_tax_gather %>%
 dplyr::group_by(condition, sgRNA) %>%
 dplyr::summarise(gene=unique(Gene),
                  mean=mean(norm_counts)) %>% # using mean expression (or median - potential outliers)
 tidyr::spread(., condition, mean)

min_expression=100
min_sgrnas=2
cat("sgRNAs before expression filtering:", nrow(norm_counts_tax_grouped), "\n") # 65959

# filtering based on expression ----
norm_counts_tax_grouped_expfilt <- norm_counts_tax_grouped %>%
  dplyr::filter( (mock >= min_expression) | (tax >= min_expression))
    
cat("sgRNAs after expression filtering:", nrow(norm_counts_tax_grouped_expfilt), "\n") # 33998

# filtering based on number of sgrnas ----
norm_counts_tax_grouped_filt <- norm_counts_tax_grouped_expfilt %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(sgrnas=paste(unique(sgRNA), collapse = ";"),
                   number_sgrnas=n()) %>%
  #dplyr::filter(number_sgrnas >= min_sgrnas) # removing control non-target
  dplyr::filter(number_sgrnas >= min_sgrnas | grepl(pattern = "NonTargetingControlGuideForMouse", gene)) # keeping control non-target

keep_sgrnas <- norm_counts_tax_grouped_filt %>%
  tidyr::separate_rows(., sgrnas, sep = ";") %>% # separate delimited values into respective rows
  dplyr::pull(sgrnas)

cat("sgRNAs after expression and sgRNA filtering:", length(keep_sgrnas), "\n")
cat("genes after expression and sgRNA filtering:", nrow(norm_counts_tax_grouped_filt), "\n")
cat("  from which control nontarget:", sum(grepl(pattern = "NonTargetingControlGuideForMouse", norm_counts_tax_grouped_filt$gene)), "\n")

# filtering raw counts and re-normalazing in mageck test!
#  on expression and number of sgRNAs and removing control non-target
raw_counts_tax_filt_exp_sgrnas_rm_nontarget <- raw_counts_tax %>%
  dplyr::filter(sgRNA %in% keep_sgrnas) %>%
  dplyr::filter(!grepl("NonTargetingControlGuideForMouse", Gene))

readr::write_tsv(raw_counts_tax_filt_exp_sgrnas_rm_nontarget, 
                 path = paste0(RESULTS_DIR, "raw_counts_tax_filt_exp_sgrnas_rm_nontarget.tsv"))