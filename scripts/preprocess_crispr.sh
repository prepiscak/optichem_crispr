#!/bin/bash

# Script to perform QC, trimming and QC after trimming
# run: preprocess_crispr.sh > preprocess_crispr.log 2>&1

NPROC=14

echo "Performing quality control step..."
# fastQC before trimming
printf '%.s*' {1..30}
printf "\nRunning fastQC before trimming: \n"
printf '%.s-' {1..20}
mkdir ./fastqc_beforeTrimming_output/
find ./ -name "*.fastq.gz" | \
    parallel -j ${NPROC} "fastqc -o ./fastqc_beforeTrimming_output/ {}"

# trimming for linked adapters
# https://github.com/marcelm/cutadapt/issues/256
# https://github.com/marcelm/cutadapt/issues/261

printf '%.s*' {1..30}
printf "\nRunning cutadapt: \n"
printf '%.s-' {1..20}
mkdir ./trimmed_fastq/

# cutadapt options: https://cutadapt.readthedocs.io/en/stable/guide.html#speed-up-tricks
# trimming linked adapters
# --discard-untrimmed - throw away untrimmed reads
# --minimum-length (-m) - discard processed reads that are shorter than LENGTH
# --maximum-length (-M) - discard processed reads that are longer than LENGTH. ?! Reads that are too long even before adapter removal are also discarded.
# for --maximum-length see also https://www.biostars.org/p/362612/

find ./ -name "*.fastq.gz" | \
    parallel -j 4 "cutadapt -j 4 -g TCTTGTGGAAAGGACGAAACACCG...GTTTTAGAGCTAGAA --discard-untrimmed -m 20 -M 20 -o {=s/_001.fastq.gz/_trimmed.fq.gz/=} {}" 
mv ./*_trimmed.fq.gz ./trimmed_fastq/

# fastQC after trimming
printf '%.s*' {1..30}
printf "\nRunning fastQC after trimming: \n"
printf '%.s-' {1..20}
mkdir ./fastqc_afterTrimming_output/
find ./ -name "*_trimmed.fq.gz" | \
    parallel -j ${NPROC} "fastqc -o ./fastqc_afterTrimming_output/ {}"                                                                   
