#!/usr/bin/bash
HUMAN_GFF_ONLY_GENES="/home/mahajrod/Reference_genomes/homo_sapiens/gff/GCF_000001405.31_GRCh38.p5_genomic_no_comments_only_chromosomes_renamed_only_genes.gff"
OUTPUT_PREFIX="/home/mahajrod/Reference_genomes/homo_sapiens/gff/GCF_000001405.31_GRCh38.p5_genomic_no_comments_only_chromosomes_renamed_only_genes"

for GENE_TYPE in antisense_RNA C_region C_region_pseudogene D_segment guide_RNA J_segment J_segment_pseudogene lncRNA miRNA misc_RNA other protein_coding pseudogene RNase_MRP_RNA RNase_P_RNA rRNA snoRNA snRNA SRP_RNA telomerase_RNA tRNA vault_RNA V_segment V_segment_pseudogene Y_RNA;
    do
    echo ${GENE_TYPE}

    grep -P "gene_biotype=${GENE_TYPE}" ${HUMAN_GFF_ONLY_GENES} > ${OUTPUT_PREFIX}.${GENE_TYPE}.gff
    done
