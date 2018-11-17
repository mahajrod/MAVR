#!/usr/bin/env bash

mkdir -p orthologs; cd orthologs; ~/Soft/eggnog-mapper-1.0.3/emapper.py -d veNOG -o amazona_agilis --report_orthologs --usemem --cpu 25 -i ../amazona_agilis.augustus.supported.transcripts.swissprot_or_pfam_and_hints.longest_pep.pep

~/Soft/eggnog-mapper-1.0.3/emapper.py -d veNOG -o amazona_collaria --report_orthologs --usemem --cpu 25 -i ../amazona_collaris.augustus.supported.transcripts.swissprot_or_pfam_and_hints.longest_pep.pep





~/Soft/MAVR/scripts/annotation/parallel_augustus.py -i JAMY3-supernova_TDTX972_v2.extended_repeat_set.selected_repeat_classes.masked.fasta -o amazona_collaria.augustus.extended_repeat_set -t 18 -x AMACOL -s chicken -c ~/Soft/augustus-3.2.1/config/ -p ~/data/db/Pfam/Pfam32.0/Pfam-A.hmm -w ~/data/db/SwissProt/SwissProt --softmasking --hintsfile all.hints.gff --extrinsicCfgFile ~/Soft/augustus-3.2.1/config/extrinsic/extrinsic.RM.EXNT.EXNS.RNASEQ.cfg


#Boqueron



for SP in cyanistes_caeruleus falco_cherrug falco_peregrinus ficedula_albicollis gallus_gallus geospiza_fortis manacus_vitellinus meleagris_gallopavo melopsittacus_undulatus parus_major serinus_canaria taeniopygia_guttata zonotrichia_albicollis; do /work/toleksyk/toleksyk/soft/MAVR/scripts/hmmer3/parallel_hmmscan.py -i /work/toleksyk/toleksyk/data/db/TreeFam/treefam9.hmm3 -s /work/toleksyk/toleksyk/longest_protein/${SP}.longest_pep.pep -o ${SP}.treefam --no_ali -t 1000 --tblout_dir splited_tblout --domtblout_dir splited_domtblout --pfamtblout_dir splited_pfamtblout -m slurm -j ${SP} -l /work/toleksyk/toleksyk/logs/${SP} -e /work/toleksyk/toleksyk/errors/${SP} -x 300 -a 120:00:00; done

melopsittacus_undulatus parus_major serinus_canaria taeniopygia_guttata zonotrichia_albicollis


#amazona_vittata annotation
cd /home/projects/carribean_parrots/amazona_vittata/skliver/genome_denovo/annotation/protein_coding_genes/fermi_SSPACE_gapcloser_LRNA_scaf/augustus_jam_repeats

cat *.hints.gff | sort -n -k 4,4 | sort -s -n -k 5,5 | sort -s -k 3,3 | sort -s -k 1,1 | join_mult_hints.pl > all.hints.gff

for I in 1 2 3; do  ~/Soft/MAVR/scripts/annotation/prepare_intron_hints_from_STAR_junction_file.py -i SJ.out.tab  -o SJ.out.m${I}.hints -m $I &  done