#!/usr/bin/env bash
#dell
BED_FILE_WITH_WINDOWS="/mnt/guatemala/skliver/1000_genomes/reference/fasta/GRCh38_full_analysis_set_plus_decoy_hla_only_chr_w_1M_s_10K.bed"
GFF_FILE_WITH_PROTEIN_CODING_GENES="/mnt/guatemala/skliver/1000_genomes/reference/gff/selected_for_GAITWAY/GCF_000001405.31_GRCh38.p5_genomic_no_comments_only_chromosomes_renamed_only_genes.protein_coding.gff"
INTERSECTION_BED_FILE="1M_10K_windows_vs_protein_coding_genes.bed"

cd /mnt/guatemala/skliver/1000_genomes/workdir/
bedtools intersect -wao \
                   -a ${BED_FILE_WITH_WINDOWS} \
                   -b ${GFF_FILE_WITH_PROTEIN_CODING_GENES} \
                    | sed 's/;.*//;s/ID=//' \
                    | awk -F'\t' '{
                                   printf "%s\t%s\t%s\t%s\n", $1,$2,$3,$12
                                   }' \
                    | awk -F'\t' 'NR==1 {
                                         printf "%s\t%s\t%s\t%s",$1,$2,$3,$4;
                                         prev_chr=$1;
                                         prev_win_start=$2
                                         };
                                  NR > 1 {
                                          if ($1 == prev_chr && $2 == prev_win_start)
                                                {
                                                 printf ",%s",$4
                                                 }
                                          else {
                                                prev_chr=$1;
                                                prev_win_start=$2;
                                                printf "\n";
                                                printf $0
                                                }
                                          }'   > ${INTERSECTION_BED_FILE}

awk -F'\t' 'NR==1 {printf "%s\t%s\t%s\t%s",$1,$2,$3,$4; prev_chr=$1; prev_win_start=$2}; NR > 1 {if ($1 == prev_chr && $2 == prev_win_start) {printf ",%s",$4} else {prev_chr=$1; prev_win_start=$2; printf "\n"; printf $0}}' 1M_10K_windows_vs_protein_coding_genes.t  > 1M_10K_windows_vs_protein_coding_genes.collapsed.bed


