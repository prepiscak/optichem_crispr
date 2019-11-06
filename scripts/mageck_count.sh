#!/bin/bash

# Script to generate count matrices
# run: mageck_count.sh > mageck_count.log 2>&1

# Using trimmed FASTQ files
sgRNA_library=./Mouse_GeCKOv2_libA_09Mar2015_ForMageck.csv
# Use only A since they share sequences with B and are therefore counted only once
control_nontarget=./Mouse_GeCKOv2_libA_nontarget_09Mar2015.csv
# Same Non-Target sgRNAs are for library A and B
FASTQ_DIR=.

# plasmid
fq_plasmid=${FASTQ_DIR}/VHR8170316_S1_R1_trimmed.fq.gz

# pretreatment
fq_pretreatment1=${FASTQ_DIR}/VHR8170314_S2_R1_trimmed.fq.gz
fq_pretreatment2=${FASTQ_DIR}/VHR8170315_S2_R1_trimmed.fq.gz
fq_pretreatment3=${FASTQ_DIR}/VHR8170316_S2_R1_trimmed.fq.gz

# mock
fq_mock1=${FASTQ_DIR}/VHR8170314_S1_R1_trimmed.fq.gz
fq_mock2=${FASTQ_DIR}/VHR8170314_S3_R1_trimmed.fq.gz
fq_mock3=${FASTQ_DIR}/VHR8170314_S4_R1_trimmed.fq.gz
fq_mock4=${FASTQ_DIR}/VHR8170314_S5_R1_trimmed.fq.gz
fq_mock5=${FASTQ_DIR}/VHR8170314_S6_R1_trimmed.fq.gz
fq_mock6=${FASTQ_DIR}/VHR8170314_S7_R1_trimmed.fq.gz
fq_mock7=${FASTQ_DIR}/VHR8170314_S8_R1_trimmed.fq.gz
fq_mock8=${FASTQ_DIR}/VHR8170314_S9_R1_trimmed.fq.gz
fq_mock9=${FASTQ_DIR}/VHR8170314_S10_R1_trimmed.fq.gz  

# tax
fq_tax1=${FASTQ_DIR}/VHR8170315_S5_R1_trimmed.fq.gz
fq_tax2=${FASTQ_DIR}/VHR8170315_S6_R1_trimmed.fq.gz
fq_tax3=${FASTQ_DIR}/VHR8170315_S7_R1_trimmed.fq.gz
fq_tax4=${FASTQ_DIR}/VHR8170316_S3_R1_trimmed.fq.gz
fq_tax5=${FASTQ_DIR}/VHR8170316_S4_R1_trimmed.fq.gz   

#mkdir ${OUTPUT_DIR}

# total normalization plasmid, pretreatment, mock, tax
# collapsing technical replicates (pre-treatment)
OUTPUT_ALL_PREFIX="total_all_tax_mock"
OUTPUT_ALL_DIR="total_all_tax_mock_counts"
SAMPLE_LABELS_ALL="plasmid_VH21,\
pretrt_VH2VH12VH22,\
mock_VH1_IL,mock_VH3_I2L2R,mock_VH4_II2R,mock_VH5_IIL,mock_VH6_IIR,mock_VH7_IIA2R,mock_VH8_IIALR,mock_VH9_IIAL,mock_VH10_I2L,\
tax_VH15_IBR,tax_VH16_IBL,tax_VH17_IB2L,tax_VH23_IBLR,tax_VH24_IB2R"

mkdir ${OUTPUT_ALL_DIR}
mageck count -l ${sgRNA_library} \
    --norm-method total \
    --output-prefix ${OUTPUT_ALL_PREFIX} \
    --sample-label ${SAMPLE_LABELS_ALL} \
    --fastq ${fq_plasmid} ${fq_pretreatment1},${fq_pretreatment2},${fq_pretreatment3} \
    ${fq_mock1} ${fq_mock2} ${fq_mock3} ${fq_mock4} ${fq_mock5} ${fq_mock6} ${fq_mock7} ${fq_mock8} ${fq_mock9} \
    ${fq_tax1} ${fq_tax2} ${fq_tax3} ${fq_tax4} ${fq_tax5}
mv ${OUTPUT_ALL_PREFIX}* ./${OUTPUT_ALL_DIR}

# only mock and tax
# total normalization mock, tax
OUTPUT_TAX_PREFIX="total_tax_mock"
OUTPUT_TAX_DIR="total_tax_mock_counts"
SAMPLE_LABELS_TAX="mock_VH1_IL,mock_VH3_I2L2R,mock_VH4_II2R,mock_VH5_IIL,mock_VH6_IIR,mock_VH7_IIA2R,mock_VH8_IIALR,mock_VH9_IIAL,mock_VH10_I2L,\
tax_VH15_IBR,tax_VH16_IBL,tax_VH17_IB2L,tax_VH23_IBLR,tax_VH24_IB2R"

mkdir ${OUTPUT_TAX_DIR}
mageck count -l ${sgRNA_library} \
    --norm-method total \
    --output-prefix ${OUTPUT_TAX_PREFIX} \
    --sample-label ${SAMPLE_LABELS_TAX} \
    --fastq ${fq_mock1} ${fq_mock2} ${fq_mock3} ${fq_mock4} ${fq_mock5} ${fq_mock6} ${fq_mock7} ${fq_mock8} ${fq_mock9} \
    ${fq_tax1} ${fq_tax2} ${fq_tax3} ${fq_tax4} ${fq_tax5}
mv ${OUTPUT_TAX_PREFIX}* ./${OUTPUT_TAX_DIR}