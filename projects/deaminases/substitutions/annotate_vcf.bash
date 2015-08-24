#!/usr/bin/env bash
WORKDIR="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/polymorphisms/"
INPUT_VCF_DIR=${WORKDIR}"vcf/"
ANNOTATED_VCF_DIR=${WORKDIR}"annotated_vcf/"
ANNOTATION_SCRIPT=~/Genetics/MAVR/scripts/annotation/annotate_vcf.py

cd ${WORKDIR}
mkdir -p ${ANNOTATED_VCF_DIR}

for SAMPLE in `ls ${INPUT_VCF_DIR}`;
    do
    SAMPLE_NAME=`basename ${SAMPLE} .vcf`
    ${ANNOTATION_SCRIPT} -i ${INPUT_VCF_DIR}${SAMPLE} -g LAN210_v0.10m -o ${ANNOTATED_VCF_DIR}${SAMPLE_NAME}_annotated.vcf \
                         -s ${ANNOTATED_VCF_DIR}${SAMPLE_NAME}.html --no_downstream --no_intron --no_upstream --no_intergenic --no_utr
    done