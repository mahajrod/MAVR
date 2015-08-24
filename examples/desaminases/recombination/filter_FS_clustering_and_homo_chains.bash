#!/usr/bin/env bash

RAW_VCF_DIR=/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/raw_vcf/
FILTER_FS_VCF_DIR=/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/raw_vcf_set_FS_filter/
FILTERED_VCF_DIR=/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/raw_vcf_strict_FS/
REFERENCE=/home/mahajrod/Genetics/Projects/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta
GATK_DIR=/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.2-0/

mkdir -p ${FILTERED_VCF_DIR} ${FILTER_FS_VCF_DIR}

cd ${RAW_VCF_DIR}

for SAMPLE in *.vcf;
    do
    java -jar ${GATK_DIR}GenomeAnalysisTK.jar -T VariantFiltration -R ${REFERENCE} -V ${SAMPLE} \
                                              --filterExpression 'FS > 20.0' --filterName 'high_strand_bias' \
                                              -o ${FILTER_FS_VCF_DIR}${SAMPLE}

    java -jar ${GATK_DIR}GenomeAnalysisTK.jar -T SelectVariants -R ${REFERENCE} -V ${FILTER_FS_VCF_DIR}${SAMPLE} \
                                              -o ${FILTERED_VCF_DIR}${SAMPLE} -ef

    done


