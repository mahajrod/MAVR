#!/usr/bin/env bash


for SP in SB6536  SB7462  SB8055  SRR1508214  SRR1508215  SRR1508749  SRR1508750; do ~/Soft/MAVR/scripts/snpcall/combine_same_sample_vcf.py -i ${SP}/splited_gvcf/splited_gvcf/ -o ${SP}/${SP}.unsorted.g.vcf -s;  ~/Soft/MAVR/scripts/snpcall/sort_vcf.py -m 200000 -p ~/Soft/picard-2.18.11/ -i ${SP}/${SP}.unsorted.g.vcf -o ${SP}/${SP}.g.vcf -e ./ -s /home/projects/mustelidae/skliver/mustela_nigripes/genome_denovo/assemblies/10x_bionano_hic/mustela_nigripes.v2.smithsonian.dict ; done