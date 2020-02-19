# Calculations of essential genes depletion
# [] - over-representation against background of all screened genes
# [] - GSEA
# [] - log2FC (median, IQR), boxplots for key genesets
# [] - Glass's delta?
# [] - use some pre-filtering for pre-treatment and plasmid, similar to what was used for treatment vs mock (100 reads)

library(devtools)
#BiocManager::install("DNAcopy")
library("DNAcopy")
#install_github("francescojm/CRISPRcleanR", force = TRUE)
library(CRISPRcleanR)
library(readr)
library(dplyr)
require("biomaRt")
require("RCurl") 
options(RCurlOptions = list(proxy="wwwcache.gla.ac.uk:8080", http.version=HTTP_VERSION_1_0)) 

# define functions ----
convertMouse2humanEntrezIDs <- function(geneList=NULL, 
                                        filters = "mgi_symbol", 
                                        ensembl_version="Ensembl Genes 96", 
                                        biomart_host="http://apr2019.archive.ensembl.org"){
  human <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl", host=biomart_host, version=ensembl_version)
  mouse <- biomaRt::useMart("ensembl", dataset="mmusculus_gene_ensembl", host=biomart_host, version=ensembl_version)
  out <- biomaRt::getLDS(attributes = c("ensembl_gene_id", "mgi_symbol"), filters = filters, values = geneList, mart = mouse, 
                         attributesL = c("ensembl_gene_id","hgnc_symbol"), 
                         martL = human,
                         verbose = TRUE)
  return(out)
}

prepareGeneFCs <- function(raw_dataset=NULL){
  
  # renaming and keeping only gene symbols, lfc, padj ----
  clean_dataset <- raw_dataset %>%
    dplyr::rename(mgi_symbol = id,
                  neg_pval = `neg|p-value`,
                  neg_padj = `neg|fdr`,
                  neg_rank = `neg|rank`,
                  neg_goodsgrna = `neg|goodsgrna`,
                  neg_lfc = `neg|lfc`) %>%
    dplyr::filter(!grepl(pattern = "NonTargetingControlGuideForMouse", mgi_symbol)) %>% # removing NonTargetingControlGuideForMouse
    dplyr::select(mgi_symbol, neg_lfc, neg_padj)
  
  #pretreatment_vs_plasmid <- pretreatment_vs_plasmid_clean %>%
  #  dplyr::select(mgi_symbol, neg_lfc, neg_padj)
  
  # Converting mouse2 human IDs ----
  dataset_human_genes_raw <- convertMouse2humanEntrezIDs(clean_dataset$mgi_symbol)
  dataset_human_genes <- dataset_human_genes_raw %>%
    dplyr::rename(mgi_symbol = MGI.symbol,
                  hgnc_symbol = HGNC.symbol) %>%
    dplyr::select(mgi_symbol, hgnc_symbol) %>%
    dplyr::filter(hgnc_symbol != "") # filtering empty human gene symbols
  
  # removing genes without human entrez_id
  clean_dataset_human <- clean_dataset %>%
    dplyr::left_join(., dataset_human_genes, by = "mgi_symbol") %>%
    dplyr::filter(!is.na(hgnc_symbol)) %>%
    dplyr::select(hgnc_symbol, neg_lfc, neg_padj) %>%
    dplyr::distinct() # removing duplicated rows
  
  # pretreatment_vs_plasmid_human_dupl <- pretreatment_vs_plasmid_human[duplicated(pretreatment_vs_plasmid_human$hgnc_symbol) |
  #                                                                       duplicated(pretreatment_vs_plasmid_human$hgnc_symbol, fromLast = TRUE),]
  
  # removing duplicates
  clean_dataset_human_LFC <- clean_dataset_human %>%
    dplyr::group_by(hgnc_symbol) %>%
    #dplyr::filter(!is.na(neg_padj)) %>% # not filtering at the moment
    # removing NA neg_padj
    dplyr::summarise(neg_lfc = neg_lfc[abs(neg_lfc) == max(abs(neg_lfc))]) %>% # keeping max absolute lfc
    dplyr::arrange(neg_lfc)
  
  return(clean_dataset_human_LFC)
}

ccr.fixedFDRthreshold_fixed <- function (FCsprofile, TruePositives, TrueNegatives, th) 
{
  presentGenes <- intersect(c(TruePositives, TrueNegatives), 
                            names(FCsprofile))
  predictions <- FCsprofile[presentGenes]
  observations <- is.element(presentGenes, TruePositives) + 
    0
  names(observations) <- presentGenes
  RES <- roc(observations, predictions, direction=">")
  COORS <- coords(RES, "all", ret = c("threshold", "ppv"), transpose = TRUE)
  FDRpercTh <- max(COORS["threshold", which(COORS["ppv", ] >= 
                                              (1 - th))])
  FDRpercRANK <- max(which(sort(FCsprofile) <= FDRpercTh))
  return(list(FCth = FDRpercTh, RANK = FDRpercRANK))
}

ccr.VisDepAndSig_fixed <- function(FCsprofile, SIGNATURES, TITLE = "", pIs = NULL, nIs = NULL, 
                                   th = 0.05, plotFCprofile = TRUE) 
{
  # fixed
  # [x] - ccr.fixedFDRthreshold_fixed; in roc() ">" changed to direction=">"
  # [x] - coords(); added transpose = TRUE after warning message for back compatibility
  # [ ] - plot not showing correct FDR based on th change!
  
  sigNames <- names(SIGNATURES)
  par(mar = c(5, 5, 8, 0))
  nsig <- length(SIGNATURES)
  if (length(pIs) > 0 & length(nIs) > 0) {
    RES <- ccr.fixedFDRthreshold_fixed(FCsprofile = FCsprofile, 
                                       TruePositives = SIGNATURES[[pIs]], TrueNegatives = SIGNATURES[[nIs]], 
                                       th = th)
    FDR5percRANK <- RES$RANK
    FDR5percTh <- RES$FCth
  } else {
    FDR5percRANK <- NULL
  }
  layout(t(matrix(c(rep(1, 5), 1:(nsig + 1), rep(1, 5), 1:(nsig + 
                                                             1)), nsig + 6, 2)))
  nelements <- length(FCsprofile)
  if (plotFCprofile) {
    plot(sort(FCsprofile, decreasing = TRUE), length(FCsprofile):1, 
         ylim = c(nelements, 1), pch = 16, frame.plot = FALSE, 
         yaxt = "n", xlim = c(min(FCsprofile), max(FCsprofile) + 
                                1), ylab = "Depletion Rank", xlab = "Log FC", 
         log = "y", cex = 1.2, cex.lab = 1.5, cex.axis = 1.2, 
         main = TITLE, cex.main = 1.5)
  }
  else {
    plot(0, 0, ylim = c(nelements, 1), pch = 16, frame.plot = FALSE, 
         yaxt = "n", xlim = c(min(FCsprofile), max(FCsprofile) + 
                                1), ylab = "Depletion Rank", xlab = "Log FC", 
         log = "y", cex = 1.2, cex.lab = 1.5, cex.axis = 1.2, 
         main = TITLE, cex.main = 1.5)
  }
  lines(x = c(0, 0), y = c(1, (length(FCsprofile)) + 10000), 
        lty = 2)
  axis(side = 2, c(1, 10, 100, 1000, 10000, nelements), labels = c(1, 
                                                                   10, 100, 1000, 10000, nelements), las = 2, cex.axis = 1)
  par(xpd = TRUE)
  lines(c(min(FCsprofile), max(FCsprofile) + 1), c(FDR5percRANK, 
                                                   FDR5percRANK), col = "red", lwd = 3, lty = 2)
  text(x = max(FCsprofile), y = FDR5percRANK, "5%FDR", pos = 3, 
       col = "red", cex = 1)
  TPR <- vector()
  N <- length(FCsprofile)
  n <- FDR5percRANK
  FTP <- vector()
  PPV <- vector()
  par(xpd = TRUE)
  for (i in 1:length(SIGNATURES)) {
    par(mar = c(5, 0, 8, 0))
    plot(1, 1, xlim = c(0, 1), ylim = c(length(FCsprofile), 
                                        1), xaxt = "n", yaxt = "n", ylab = "", xlab = "", 
         col = NA, log = "y", frame.plot = FALSE)
    hitPositions <- match(SIGNATURES[[i]], names(sort(FCsprofile)))
    hitPositions <- hitPositions[!is.na(hitPositions)]
    abline(h = hitPositions[hitPositions <= FDR5percRANK], 
           lwd = 2, col = "blue")
    abline(h = hitPositions[hitPositions > FDR5percRANK], 
           lwd = 2, col = "gray")
    TPR[i] <- length(which(hitPositions <= FDR5percRANK))/length(hitPositions)
    lines(c(0, 1), c(FDR5percRANK, FDR5percRANK), col = "red", 
          lwd = 3, lty = 2)
    text(0.4, 1, labels = sigNames[i], pos = 4, offset = 0, 
         srt = 80, cex = 1)
    if (i < (length(SIGNATURES) - 1)) {
      text(0.6, nelements + 50000, paste(round(100 * TPR[i]), 
                                         "%", sep = ""), cex = 1.2, col = "blue")
    }
  }
  names(TPR) <- sigNames
  return(TPR)
}


# Loading results of Mageck test ---- 
# geneFCs - The function ccr.geneMeanFCs performs this conversion by considering for each gene the average logFC across its targeting sgRNAs.

# preatreatment_vs_plasmid
pretreatment_vs_plasmid_raw <- readr::read_tsv("./pretreatment_vs_plasmid/pretreatment_vs_plasmid.gene_summary.txt",
                                               col_types = cols(
                                                 id = col_character(),
                                                 num = col_double(),
                                                 `neg|score` = col_double(),
                                                 `neg|p-value` = col_double(),
                                                 `neg|fdr` = col_double(),
                                                 `neg|rank` = col_double(),
                                                 `neg|goodsgrna` = col_double(),
                                                 `neg|lfc` = col_double(),
                                                 `pos|score` = col_double(),
                                                 `pos|p-value` = col_double(),
                                                 `pos|fdr` = col_double(),
                                                 `pos|rank` = col_double(),
                                                 `pos|goodsgrna` = col_double(),
                                                 `pos|lfc` = col_double()
                                               ))

pretreatment_vs_plasmid_human_LFC <- prepareGeneFCs(raw_dataset = pretreatment_vs_plasmid_raw)

geneFCs <- pretreatment_vs_plasmid_human_LFC$neg_lfc
names(geneFCs) <- pretreatment_vs_plasmid_human_LFC$hgnc_symbol

# mock_vs_pretreatment
mock_vs_pretreatment_raw <- readr::read_tsv("./mock_vs_pretreatment/mock_vs_pretreatment.gene_summary.txt", 
                                            col_types = cols(
                                              id = col_character(),
                                              num = col_double(),
                                              `neg|score` = col_double(),
                                              `neg|p-value` = col_double(),
                                              `neg|fdr` = col_double(),
                                              `neg|rank` = col_double(),
                                              `neg|goodsgrna` = col_double(),
                                              `neg|lfc` = col_double(),
                                              `pos|score` = col_double(),
                                              `pos|p-value` = col_double(),
                                              `pos|fdr` = col_double(),
                                              `pos|rank` = col_double(),
                                              `pos|goodsgrna` = col_double(),
                                              `pos|lfc` = col_double()
                                            ))

mock_vs_pretreatment_human_LFC <- prepareGeneFCs(raw_dataset = mock_vs_pretreatment_raw)

plot(mock_vs_pretreatment_human_LFC$neg_lfc, pretreatment_vs_plasmid_human_LFC$neg_lfc)
cor.test(mock_vs_pretreatment_human_LFC$neg_lfc, pretreatment_vs_plasmid_human_LFC$neg_lfc, method = "spearman")


geneFCs <- mock_vs_pretreatment_human_LFC$neg_lfc
names(geneFCs) <- mock_vs_pretreatment_human_LFC$hgnc_symbol

# tax_vs_pretreatment
tax_vs_pretreatment_raw <- readr::read_tsv("./tax_vs_pretreatment/tax_vs_pretreatment.gene_summary.txt", 
                                            col_types = cols(
                                              id = col_character(),
                                              num = col_double(),
                                              `neg|score` = col_double(),
                                              `neg|p-value` = col_double(),
                                              `neg|fdr` = col_double(),
                                              `neg|rank` = col_double(),
                                              `neg|goodsgrna` = col_double(),
                                              `neg|lfc` = col_double(),
                                              `pos|score` = col_double(),
                                              `pos|p-value` = col_double(),
                                              `pos|fdr` = col_double(),
                                              `pos|rank` = col_double(),
                                              `pos|goodsgrna` = col_double(),
                                              `pos|lfc` = col_double()
                                            ))

tax_vs_pretreatment_human_LFC <- prepareGeneFCs(raw_dataset = tax_vs_pretreatment_raw)

plot(mock_vs_pretreatment_human_LFC$neg_lfc, pretreatment_vs_plasmid_human_LFC$neg_lfc)
plot(mock_vs_pretreatment_human_LFC$neg_lfc, tax_vs_pretreatment_human_LFC$neg_lfc)
cor.test(mock_vs_pretreatment_human_LFC$neg_lfc, pretreatment_vs_plasmid_human_LFC$neg_lfc, method = "pearson")
cor.test(mock_vs_pretreatment_human_LFC$neg_lfc, tax_vs_pretreatment_human_LFC$neg_lfc, method = "pearson")
cor.test(tax_vs_pretreatment_human_LFC$neg_lfc, pretreatment_vs_plasmid_human_LFC$neg_lfc, method = "pearson")

geneFCs <- tax_vs_pretreatment_human_LFC$neg_lfc
names(geneFCs) <- tax_vs_pretreatment_human_LFC$hgnc_symbol

# when enrichment is done it seems that Ribosomal_Proteins and CFE are down! but not in mock vs vs pretreatment ( differences between mocks!?)

# treatment_vs_pretreatment
treatment_vs_pretreatment_raw <- readr::read_tsv("./treatment_vs_pretreatment/treatment_vs_pretreatment.gene_summary.txt", 
                                           col_types = cols(
                                             id = col_character(),
                                             num = col_double(),
                                             `neg|score` = col_double(),
                                             `neg|p-value` = col_double(),
                                             `neg|fdr` = col_double(),
                                             `neg|rank` = col_double(),
                                             `neg|goodsgrna` = col_double(),
                                             `neg|lfc` = col_double(),
                                             `pos|score` = col_double(),
                                             `pos|p-value` = col_double(),
                                             `pos|fdr` = col_double(),
                                             `pos|rank` = col_double(),
                                             `pos|goodsgrna` = col_double(),
                                             `pos|lfc` = col_double()
                                           ))

treatment_vs_pretreatment_human_LFC <- prepareGeneFCs(raw_dataset = treatment_vs_pretreatment_raw)

geneFCs <- treatment_vs_pretreatment_human_LFC$neg_lfc
names(geneFCs) <- treatment_vs_pretreatment_human_LFC$hgnc_symbol


# tax_vs_mock
tax_vs_mock_raw <- readr::read_tsv("/media/prepisca/DataAnalysis1/optichem_crispr/optichem_crispr/datasets/tax_mock_test/tax_mock_test.gene_summary.txt", 
                                           col_types = cols(
                                             id = col_character(),
                                             num = col_double(),
                                             `neg|score` = col_double(),
                                             `neg|p-value` = col_double(),
                                             `neg|fdr` = col_double(),
                                             `neg|rank` = col_double(),
                                             `neg|goodsgrna` = col_double(),
                                             `neg|lfc` = col_double(),
                                             `pos|score` = col_double(),
                                             `pos|p-value` = col_double(),
                                             `pos|fdr` = col_double(),
                                             `pos|rank` = col_double(),
                                             `pos|goodsgrna` = col_double(),
                                             `pos|lfc` = col_double()
                                           ))

tax_vs_mock_human_LFC <- prepareGeneFCs(raw_dataset = tax_vs_mock_raw)

tax_vs_mock_geneFCs <- tax_vs_mock_human_LFC$neg_lfc
names(tax_vs_mock_geneFCs) <- tax_vs_mock_human_LFC$hgnc_symbol

tax_vs_mock_recall_scores <- ccr.VisDepAndSig_fixed(FCsprofile = tax_vs_mock_geneFCs,
                                        SIGNATURES = SIGNATURES,
                                        TITLE = '',
                                        pIs = 6,
                                        nIs = 7,
                                        th = 0.5,
                                        plotFCprofile=TRUE)

tax_vs_mock_recall_scores

# enrichment testing ----
data(EssGenes.ribosomalProteins)
data(EssGenes.DNA_REPLICATION_cons)
data(EssGenes.KEGG_rna_polymerase)
data(EssGenes.PROTEASOME_cons)
data(EssGenes.SPLICEOSOME_cons)
data(BAGEL_essential)
data(BAGEL_nonEssential)

SIGNATURES <- list(Ribosomal_Proteins=EssGenes.ribosomalProteins, 
                   DNA_Replication = EssGenes.DNA_REPLICATION_cons, 
                   RNA_polymerase = EssGenes.KEGG_rna_polymerase, 
                   Proteasome = EssGenes.PROTEASOME_cons, 
                   Spliceosome = EssGenes.SPLICEOSOME_cons, 
                   CFE = BAGEL_essential, 
                   non_essential = BAGEL_nonEssential)

Ribosomal_Proteins_df <- data.frame(geneset = "Ribosomal_Proteins",
                                    genes = EssGenes.ribosomalProteins)
DNA_Replication_df <- data.frame(geneset = "DNA_Replication",
                                    genes = EssGenes.DNA_REPLICATION_cons)
RNA_polymerase_df <- data.frame(geneset = "RNA_polymerase",
                                    genes = EssGenes.KEGG_rna_polymerase)
Proteasome_df <- data.frame(geneset = "Proteasome",
                                genes = EssGenes.PROTEASOME_cons)
Spliceosome_df <- data.frame(geneset = "Spliceosome",
                            genes = EssGenes.SPLICEOSOME_cons)
essential_df <- data.frame(geneset = "essential",
                             genes = BAGEL_essential)
non_essential_df <- data.frame(geneset = "non_essential",
                     genes = BAGEL_nonEssential)

SIGNATURES_df <- bind_rows(Ribosomal_Proteins_df, DNA_Replication_df, RNA_polymerase_df,
                           Proteasome_df, Spliceosome_df, essential_df, non_essential_df) %>%
  dplyr::rename(hgnc_symbol = genes)


pretreatment_vs_plasmid_human_signatures <- mock_vs_pretreatment_human_LFC %>%
  dplyr::left_join(., SIGNATURES_df, by = "hgnc_symbol") %>%
  dplyr::filter(!is.na(geneset))
pretreatment_vs_plasmid_human_signatures <- pretreatment_vs_plasmid_human_signatures %>%
  dplyr::mutate(geneset = factor(geneset, 
                                 levels = c("non_essential", "essential",
                                            "Ribosomal_Proteins", "Proteasome", "Spliceosome",
                                            "DNA_Replication", "RNA_polymerase"))) %>%
  dplyr::rename(log2FC = neg_lfc)
  
  
library(ggplot2)

pretreatment_vs_plasmid_human_signatures %>%
  group_by(geneset) %>%
  arrange(., .by_group = TRUE) %>%
  ggplot(.) +
  geom_point(aes(x=hgnc_symbol, y=log2FC)) +
  facet_wrap(~geneset, scales = "free") +
  theme_minimal()

pretreatment_vs_plasmid_human_signatures %>%
  dplyr::mutate(geneset = factor(geneset, 
                                 levels = c("non_essential", "essential",
                                            "Ribosomal_Proteins", "Proteasome", "Spliceosome",
                                            "DNA_Replication", "RNA_polymerase"))) 

ggpubr::ggviolin(pretreatment_vs_plasmid_human_signatures,
                  "geneset", "log2FC", 
                  fill = "geneset", 
                 add = "boxplot", add.params = list(fill = "white")) + 
  geom_hline(yintercept = 0)

# test 
#data(EPLC.272HcorrectedFCs)
# ## storing the sgRNA log fold changes into a name vector
#FCs<-EPLC.272HcorrectedFCs$corrected_logFCs$avgFC
#names(FCs)<-rownames(EPLC.272HcorrectedFCs$corrected_logFCs)
# ## loading sgRNA library annotation
data(KY_Library_v1.0)
# ## computing gene average log fold changes
#FCs<-ccr.geneMeanFCs(FCs,KY_Library_v1.0)

# manual calculation of mean gene FC
pretreatment_vs_plasmid_sgrnas_raw <- readr::read_tsv("./pretreatment_vs_plasmid/pretreatment_vs_plasmid.sgrna_summary.txt")

pretreatment_vs_plasmid_sgrnas <- pretreatment_vs_plasmid_sgrnas_raw %>%
  dplyr::filter(!grepl(pattern = "NonTargetingControlGuideForMouse", Gene)) # removing NonTargetingControlGuideForMouse

pretreatment_vs_plasmid_sgrnas_library <- pretreatment_vs_plasmid_sgrnas %>%
  dplyr::select(sgrna, Gene) %>%
  tibble::column_to_rownames(., var="sgrna")

pretreatment_vs_plasmid_sgrnas_LFC <- pretreatment_vs_plasmid_sgrnas$LFC
names(pretreatment_vs_plasmid_sgrnas_LFC) <- pretreatment_vs_plasmid_sgrnas$sgrna

geneMeanFCs <- function(sgRNA_FCprofile, libraryAnnotation) 
{
  FCsprofile <- sgRNA_FCprofile
  FCsprofile <- aggregate(sgRNA_FCprofile ~ libraryAnnotation[names(sgRNA_FCprofile), 
                                                              "Gene"], FUN = mean)
  nn <- as.character(FCsprofile$`libraryAnnotation[names(sgRNA_FCprofile), "Gene"]`)
  FCsprofile <- FCsprofile$sgRNA_FCprofile
  names(FCsprofile) <- as.character(nn)
  return(FCsprofile)
}

geneFCs2 <- geneMeanFCs(pretreatment_vs_plasmid_sgrnas_LFC, pretreatment_vs_plasmid_sgrnas_library)

# Converting mouse2 human IDs ----
dataset_human_genes_raw <- convertMouse2humanEntrezIDs(pretreatment_vs_plasmid_sgrnas$Gene)
dataset_human_genes <- dataset_human_genes_raw %>%
  dplyr::rename(mgi_symbol = MGI.symbol,
                hgnc_symbol = HGNC.symbol) %>%
  dplyr::select(mgi_symbol, hgnc_symbol) %>%
  dplyr::filter(hgnc_symbol != "") # filtering empty human gene symbols

# removing genes without human entrez_id
clean_dataset_human <- pretreatment_vs_plasmid_sgrnas %>%
  dplyr::left_join(., dataset_human_genes, by = c("Gene"="mgi_symbol")) %>%
  dplyr::filter(!is.na(hgnc_symbol)) %>%
  dplyr::select(sgrna, hgnc_symbol, LFC) %>%
  dplyr::distinct() # removing duplicated rows

# pretreatment_vs_plasmid_human_dupl <- pretreatment_vs_plasmid_human[duplicated(pretreatment_vs_plasmid_human$hgnc_symbol) |
#                                                                       duplicated(pretreatment_vs_plasmid_human$hgnc_symbol, fromLast = TRUE),]

# removing duplicates
clean_dataset_human_LFC <- clean_dataset_human %>%
  dplyr::group_by(hgnc_symbol) %>%
  #dplyr::filter(!is.na(neg_padj)) %>% # not filtering at the moment
  # removing NA neg_padj
  dplyr::summarise(neg_lfc = neg_lfc[abs(neg_lfc) == max(abs(neg_lfc))]) %>% # keeping max absolute lfc
  dplyr::arrange(neg_lfc)

geneFCs2_df <- data.frame(mgi_symbol = names(geneFCs2), LFC=geneFCs2) %>%
  dplyr::left_join(., dataset_human_genes, by = "mgi_symbol") %>%
  dplyr::select(-mgi_symbol) %>%
  dplyr::filter(!is.na(hgnc_symbol)) %>%
  dplyr::group_by(hgnc_symbol) %>%
  dplyr::summarise(LFC = unique(LFC[abs(LFC) == max(abs(LFC))])) %>% # keeping max absolute lfc
  dplyr::arrange(LFC)
  
#geneFCs2_df_dupl <- geneFCs2_df[duplicated(geneFCs2_df$hgnc_symbol) | duplicated(geneFCs2_df$hgnc_symbol, fromLast = TRUE),]


geneFCs2 <- geneFCs2_df$LFC
names(geneFCs2) <- geneFCs2_df$hgnc_symbol
  
## Visualising log fold change profile with superimposed signatures specifying
## that the reference gene sets are in positions 6 and 7
Recall_scores <- ccr.VisDepAndSig_fixed(FCsprofile = geneFCs,
                                SIGNATURES = SIGNATURES,
                                TITLE = '',
                                pIs = 6,
                                nIs = 7,
                                th = 0.7,
                                plotFCprofile=TRUE)

Recall_scores
