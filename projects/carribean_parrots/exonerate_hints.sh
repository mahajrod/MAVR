#!/usr/bin/env bash

cd /home/projects/carribean_parrots/amazona_collaria/skliver/genome_denovo/annotation/protein_coding_genes/S2/exonerate/

~/Soft/MAVR/scripts/annotation/split_exonerate_output.py -i gallus_gallus/out -g AMACOL.GALGAL.G -o gallus_gallus/amazona_collaria.gallus_gallus -t AMACOL.GALGAL.T -r ~/data/db/proteins/gallus_gallus/ncbi/gallus_gallus.no_custom_aa.pep &
~/Soft/MAVR/scripts/annotation/split_exonerate_output.py -i melopsittacus_undulatus/out -g AMACOL.MELUND.G -o melopsittacus_undulatus/amazona_collaria.melopsittacus_undulatus -t AMACOL.MELUND.T -r ~/data/db/proteins/melopsittacus_undulatus/ncbi/melopsittacus_undulatus.no_custom_aa.pep &
~/Soft/MAVR/scripts/annotation/split_exonerate_output.py -i taeniopygia_guttata/out -g AMACOL.TAEGUT.G -o taeniopygia_guttata/amazona_collaria.taeniopygia_guttata -t AMACOL.TAEGUT.T -r ~/data/db/proteins/taeniopygia_guttata/ncbi/taeniopygia_guttata.no_custom_aa.pep &

for SP in gallus_gallus taeniopygia_guttata melopsittacus_undulatus;
    do
    ~/Soft/MAVR/scripts/annotation/prepare_hints_from_exonerate_target_output.py -i ${SP}/amazona_collaria.${SP}.target.gff  \
                                                                                 -o ${SP}/amazona_collaria.${SP}.target \
                                                                                 -m 5 \
                                                                                 -e ~/Soft/augustus-3.2.1/scripts/ &
    done



~/Soft/MAVR/scripts/annotation/split_exonerate_output.py -i gallus_gallus/out -g AMAAGI.GALGAL.G -o gallus_gallus/amazona_agilis.gallus_gallus -t AMAAGI.GALGAL.T -r ~/data/db/proteins/gallus_gallus/ncbi/gallus_gallus.no_custom_aa.pep &
~/Soft/MAVR/scripts/annotation/split_exonerate_output.py -i melopsittacus_undulatus/out -g AMAAGI.MELUND.G -o melopsittacus_undulatus/amazona_agilis.melopsittacus_undulatus -t AMAAGI.MELUND.T -r ~/data/db/proteins/melopsittacus_undulatus/ncbi/melopsittacus_undulatus.no_custom_aa.pep &
~/Soft/MAVR/scripts/annotation/split_exonerate_output.py -i taeniopygia_guttata/out -g AMAAGI.TAEGUT.G -o taeniopygia_guttata/amazona_agilis.taeniopygia_guttata -t AMAAGI.TAEGUT.T -r ~/data/db/proteins/taeniopygia_guttata/ncbi/taeniopygia_guttata.no_custom_aa.pep &

for SP in gallus_gallus taeniopygia_guttata melopsittacus_undulatus;
    do
    ~/Soft/MAVR/scripts/annotation/prepare_hints_from_exonerate_target_output.py -i ${SP}/amazona_agilis.${SP}.target.gff  \
                                                                                 -o ${SP}/amazona_agilis.${SP}.target \
                                                                                 -m 5 \
                                                                                 -e ~/Soft/augustus-3.2.1/scripts/ &
    done
