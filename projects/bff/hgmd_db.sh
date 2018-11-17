#!/usr/bin/env bash

cd /home/projects/mustelidae/skliver/analysis/HGMD/bff_snps_in_hgmd_genes/

for ID in `ls ../hgmd_variants/ | grep gene.ids`; do cat /home/projects/mustelidae/skliver/mustela_nigripes/genome_ref_assisted/SNPcall/variant_annotation/7sp_merged.sorted.snpeff.combined.good.no_RM_repeats.common_names.snpeff | grep -v downstream_gene_variant | grep -v upstream_gene_variant| grep -v intergenic_region | grep -v intron_variant | ~/Soft/MAVR/scripts/file/extract_by_column_value.py -f ../hgmd_variants/${ID} -c "-1" -o gene_variants_no_intron/7sp_merged.sorted.snpeff.combined.good.no_RM_repeats.common_names.${ID%.ids}.no_intron.snpeff & done

for FILE in `ls gene_variants_no_intron`; do awk -F'\t' 'NR == 1 {print $0}; NR > 1 {printf "%s\t%s\t%s\t%s\n",$22,$15,$16,$7}' gene_variants_no_intron/${FILE} > gene_variants_no_intron_short/${FILE%snpeff}short.snpeff & done



cd /home/projects/mustelidae/skliver/analysis/HGMD/

for FILE in `ls hgmd_variants| grep txt`; do awk -F'\t' 'NR == 1 {printf "#%s\t%s\n",$4,$8}; NR > 1 {printf "%s\t%s\n",$4,$8}' hgmd_variants/${FILE} > hgmd_variants_short/${FILE%txt}short  & done


mkdir hgmd_variants_short_collapsed/