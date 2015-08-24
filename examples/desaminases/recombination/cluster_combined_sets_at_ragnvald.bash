#!/usr/bin/env bash

MACE_DIR="/home/mahajrod/Genetics/MACE/"
SCRIPT=${MACE_DIR}/scripts/clustering_pipeline.py

WORKDIR=/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/cluster_raw_vcf_strict_FS/

VCF_DIR="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/cluster_raw_vcf_strict_FS/combined/"
REFERENCE=/home/mahajrod/Genetics/Projects/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta
REFERENCE_ANNOTATIONS=/home/mahajrod/Genetics/Projects/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3
REFERENCE_MASKING=/home/mahajrod/Genetics/Projects/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all_not_in_good_genes.gff

SAMPLE_SETS=(A1_D1 A1_D3 A1_D6 AID_D1 AID_D3 AID_D6 PmCDA1_D1 PmCDA1_D3 PmCDA1_D6)

cd ${WORKDIR}

for SET in ${SAMPLE_SETS[@]};
    do
    CURRENT_SAMPLE_SET=${SET}
    ${SCRIPT} -r ${REFERENCE} -a ${REFERENCE_ANNOTATIONS} -m ${REFERENCE_MASKING} -f ${VCF_DIR}${CURRENT_SAMPLE_SET}_raw.vcf -s ${CURRENT_SAMPLE_SET} -y ${CURRENT_SAMPLE_SET}
    done