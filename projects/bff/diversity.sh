#!/usr/bin/env bash

cd /mnt/47/projects/mustelidae/skliver/mustela_nigripes/genome_ref_assisted/SNPcall/10x_bionano_hic/vcf/all/

 ~/Soft/MACE/scripts/pandas/draw_zygoty.py -i bff_7_samples.snp.good.no_repeats.vcf -o bff_7_samples.snp.good.no_repeats.zygoty

 ~/Soft/MACE/scripts/pandas/draw_zygoty.py -i bff_7_samples.indel.good.no_repeats.vcf -o bff_7_samples.indel.good.no_repeats.zygoty
