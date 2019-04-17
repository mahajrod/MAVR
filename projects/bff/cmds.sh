#!/usr/bin/env bash


for SP in SB6536; do mkdir -p ${SP}/splited_gvcf; cd ${SP}/splited_gvcf; ~/Soft/MAVR/scripts/snpcall/parallel_gvcf_call.py -b /home/projects/mustelidae/skliver/mustela_nigripes/genome_ref_assisted/alignment/10x_bionano_hic/${SP}/${SP}.mkdup.bam -o ./ -p ${SP} -r /home/projects/mustelidae/skliver/mustela_nigripes/genome_denovo/assemblies/10x_bionano_hic/mustela_nigripes.v2.smithsonian.fasta -g ~/Soft/GenomeAnalysisTK-3.7/ -t 40 -m 280g -x 15 -l 300000 ; cd ../../; done


for SP in SB6536  SB7462  SB8055  SRR1508214  SRR1508215  SRR1508749  SRR1508750; do ~/Soft/MAVR/scripts/snpcall/combine_same_sample_vcf.py -i ${SP}/splited_gvcf/splited_gvcf/ -o ${SP}/${SP}.unsorted.g.vcf -s;  ~/Soft/MAVR/scripts/snpcall/sort_vcf.py -m 200000 -p ~/Soft/picard-2.18.11/ -i ${SP}/${SP}.unsorted.g.vcf -o ${SP}/${SP}.g.vcf -e ./ -s /home/projects/mustelidae/skliver/mustela_nigripes/genome_denovo/assemblies/10x_bionano_hic/mustela_nigripes.v2.smithsonian.dict ; done


for sample in `bcftools query -l ../../all/bff_7_samples.snp.good.vcf`; do java -jar ~/Soft/GenomeAnalysisTK-3.7/GenomeAnalysisTK.jar -T SelectVariants -R /home/projects/mustelidae/skliver/mustela_nigripes//genome_denovo/assemblies/10x_bionano_hic/mustela_nigripes.v2.smithsonian.fasta -V ../../all/bff_7_samples.snp.good.vcf -o ${sample}.snp.good.vcf -sn $sample -env -ef ; done

~/Soft/MACE/scripts/draw_variant_window_densities.py -i ../../SB6536.snp.good.hetero.vcff -r /home/projects/mustelidae/skliver/mustela_nigripes/genome_denovo/assemblies/10x_bionano_hic/mustela_nigripes.v2.smithsonian.fasta -p parse -a HiC_scaffold_15,HiC_scaffold_12,HiC_scaffold_13,HiC_scaffold_14,HiC_scaffold_5,HiC_scaffold_18,HiC_scaffold_10,HiC_scaffold_6,HiC_scaffold_8,HiC_scaffold_11,HiC_scaffold_7,HiC_scaffold_17,HiC_scaffold_2,HiC_scaffold_9,HiC_scaffold_19,HiC_scaffold_1,HiC_scaffold_3,HiC_scaffold_16,HiC_scaffold_4 -z HiC_scaffold_15,HiC_scaffold_12,HiC_scaffold_13,HiC_scaffold_14,HiC_scaffold_5,HiC_scaffold_18,HiC_scaffold_10,HiC_scaffold_6,HiC_scaffold_8,HiC_scaffold_11,HiC_scaffold_7,HiC_scaffold_17,HiC_scaffold_2,HiC_scaffold_9,HiC_scaffold_19,HiC_scaffold_1,HiC_scaffold_3,HiC_scaffold_16,HiC_scaffold_4 -o SB6536.100k_jet.10cat --masking_threshold 0.5 -w 100000 --colormap jet --density_thresholds 0.0,0.05,0.1,0.25,0.5,0.75,1.0,1.5,2.0,2.5




for SP in enhydra_lutris mustela_nigripes neovison_vison homo_sapiens; do   mkdir -p predicted_gene_names; awk -F'\t' -v SP="${SP}" '{printf "%s\t%s@%s\n", $2,SP,  $1}' ${SP}.predicted_gene_names > predicted_gene_names/${SP}.predicted_gene_names & done
mkdir predicted_gene_names; mv *.predicted_gene_names predicted_gene_names/
for SP in enhydra_lutris mustela_nigripes neovison_vison homo_sapiens; do ~/Soft/MAVR/scripts/sequence/label_sequences.py -i pep/${SP}.pep -o labeled_pep/${SP}.pep -l ${SP};  done
for SP in enhydra_lutris mustela_nigripes neovison_vison homo_sapiens; do sed "s/^/${SP}@/;s/\t/\t${SP}@/" accordance/${SP}.cds_to_pep.accordance > labeled_accordance/${SP}.cds_to_pep.accordance & done
for SP in enhydra_lutris mustela_nigripes neovison_vison homo_sapiens; do ~/Soft/MAVR/scripts/sequence/label_sequences.py -i cds/${SP}.cds -o labeled_cds/${SP}.cds -l ${SP};  done