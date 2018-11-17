#!/usr/bin/env bash


for SP in amazona_leucocephala amazona_ventralis amazona_vittata aquila_chrysaetos cyanistes_caeruleus falco_cherrug falco_peregrinus ficedula_albicollis geospiza_fortis manacus_vitellinus melopsittacus_undulatus parus_major serinus_canaria taeniopygia_guttata zonotrichia_albicollis; do ~/Soft/MAVR/scripts/annotation/trim_cds_and_remove_terminal_stop_codons.py -i cds/${SP}.with_pep.cds -o trimmed_cds/${SP}.trimmed.cds ; done



for SET in set2 set3 set4 set5; do ~/soft/MAVR/scripts/evolution/divergence_time_estimation_all_clock.py -s ${SET}.4fold_degenerated_sites.fasta -t ${SET}.nwk -o ${SET}_200k_it2 -m HKY85 -d ~/soft/PAML/paml4.7a/src/ --number_of_samples 200000 --num_of_burning 20000 &  done


for SP in amazona_leucocephala amazona_ventralis amazona_vittata; do bedtools intersect -v -header -a ${SP}.snp.vcf -b ${SP}.uncallable.gff > ${SP}.snp.callable_only.vcf &  done


