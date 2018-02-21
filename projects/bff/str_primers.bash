#!/usr/bin/env bash

cd ~/projects/mustelidae/mustela_nigripes/genome_denovo/annotation/repeats/hybrid_assembly/TRF

~/Soft/MAVR/scripts/repeat_masking/tandem_repeat_masking.py -i ../../../../assemblies/bionano/assemblies/hybrid_assembly/assembly.hybrid.all.fasta -o assembly.hybrid.all -t 30 -p ~/Soft/TRF/trf

~/Soft/MAVR/scripts/repeat_masking/filter_trf_gff.py -i assembly.hybrid.all.trf.with_rep_seqs.gff -o assembly.hybrid.all.trf.with_rep_seqs.monomer_4.copy_no_less_20.gff -x assembly.hybrid.all.trf.with_rep_seqs.monomer_4.copy_no_less_20.filtered_out.gff -n 4 -m 4 -b 20

~/Soft/MAVR/scripts/repeat_masking/filter_trf_gff_by_exact_copy_number.py -i assembly.hybrid.all.trf.with_rep_seqs.monomer_4.copy_no_less_20.gff -o assembly.hybrid.all.trf.with_rep_seqs.monomer_4.exact_tandem_copy_no_less_20.gff -x assembly.hybrid.all.trf.with_rep_seqs.monomer_4.exact_tandem_copy_no_less_20.filtered_out.gff -b 20 -p

~/Soft/MAVR/scripts/repeat_masking/filter_trf_gff_by_exact_copy_number.py -i assembly.hybrid.all.trf.with_rep_seqs.monomer_4.copy_no_less_20.gff -o assembly.hybrid.all.trf.with_rep_seqs.monomer_4.exact_copy_no_less_20.gff -x assembly.hybrid.all.trf.with_rep_seqs.monomer_4.exact_copy_no_less_20.filtered_out.gff -b 20

~/Soft/MAVR/scripts/annotation/gff/add_flanks_to_gff_record.py  -i assembly.hybrid.all.trf.with_rep_seqs.monomer_4.exact_tandem_copy_no_less_20.gff -o assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks -f ../../../../assemblies/bionano/assemblies/hybrid_assembly/assembly.hybrid.all.fasta

~/Soft/MAVR/scripts/sequence/extract_sequences_by_gff.py -i ../../../../assemblies/bionano/assemblies/hybrid_assembly/assembly.hybrid.all.fasta -g assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.gff -o assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.fasta -t repeat

~/Soft/MAVR/scripts/primers/predict_str_primers.py -f assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.gff -s assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.fasta -o assembly.hybrid.all.trf -k ../../../../assemblies/bionano/assemblies/hybrid_assembly/kmers/ -r mustela_nigripes_hybrid -p ~/Soft/primer3/src/ -y ~/Soft/primer3/src/primer3_config/ -m

~/Soft/MAVR/scripts/primers/predict_str_primers.py -f assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.gff -s assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.fasta -o assembly.hybrid.all.trf -k ../../../../assemblies/bionano/assemblies/hybrid_assembly/kmers/ -r mustela_nigripes_hybrid -p ~/Soft/primer3/src/ -y ~/Soft/primer3/src/primer3_config/
