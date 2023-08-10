#!/usr/bin/env bash

conda install -y  mamba
mamba install -y biopython numpy scipy pandas matplotlib venn xmltodict bcbio-gff statsmodels xlsxwriter ete3

mamba install -y fastqc trimmomatic bwa samtools bowtie2 blast mosdepth