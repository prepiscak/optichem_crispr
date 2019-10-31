#!/bin/bash

# Script to perform mageck test on raw pre-filtered count tables

# tax
TREATMENT_IDs="tax_VH15_IBR,tax_VH16_IBL,tax_VH17_IB2L,tax_VH23_IBLR,tax_VH24_IB2R"
# mock
CONTROL_IDs="mock_VH1_IL,mock_VH3_I2L2R,mock_VH4_II2R,mock_VH5_IIL,mock_VH6_IIR,mock_VH7_IIA2R,mock_VH8_IIALR,mock_VH9_IIAL,mock_VH10_I2L"

# total normalization mock, tax - filtering expression and number of sgRNAs and removing non-target
# using raw counts after filtering and re-normalizing
OUTPUT_PREFIX="tax_mock_renorm"
OUTPUT_DIR="tax_mock_renorm"
count_table="./crispr_screen_filter/raw_counts_tax_filt_exp_sgrnas_rm_nontarget.tsv"

mkdir ${OUTPUT_DIR}
mageck test \
    --norm-method total \
    --adjust-method fdr \
    --normcounts-to-file --keep-tmp \
    --count-table ${count_table} \
    --treatment-id ${TREATMENT_IDs}\
    --control-id ${CONTROL_IDs}\
    --output-prefix ${OUTPUT_PREFIX} \
    --additional-rra-parameters "--min-percentage-goodsgrna 0.6" 

mv ${OUTPUT_PREFIX}* ./${OUTPUT_DIR}
