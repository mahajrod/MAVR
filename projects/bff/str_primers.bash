#!/usr/bin/env bash

cd ~/projects/mustelidae/mustela_nigripes/genome_denovo/annotation/repeats/hybrid_assembly/TRF

~/Soft/MAVR/scripts/repeat_masking/tandem_repeat_masking.py -i ../../../../assemblies/bionano/assemblies/hybrid_assembly/assembly.hybrid.all.fasta -o assembly.hybrid.all -t 30 -p ~/Soft/TRF/trf


