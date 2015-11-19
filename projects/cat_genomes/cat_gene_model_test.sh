#!/usr/bin/env bash
#supermicro
cd workdir/cat/cat

sed -r 's/.*Class=(.*);Family.*/\1/' scaf_fca70.fa.gff | sort | uniq > scaf_fca70.fa.annotated_repeat_types.t
grep -P "Class=DNA|Class=DNA\?|Class=LINE|Class=LTR|Class=LTR\?|Class=RC|Class=RC\?|Class=SINE|Class=SINE\?" scaf_fca70.fa.gff > scaf_fca70.fa.selected_repeat_types.gff

grep -P "\tCDS\t" augustus.gff > augustus_only_CDS.gff

bedtools intersect -u -a augustus_only_CDS.gff \
                   -b ~/data/genomes/cat/masking/scaf_fca70.fa.selected_repeat_types.gff > augustus_only_CDS.gff.intersecting_with_repeats.gff

sed "s/.*=//" augustus_only_CDS.gff.intersecting_with_repeats.gff | sort | uniq > augustus_only_CDS.gff.intersecting_with_repeats.ids

~/soft/MAVR/scripts/annotation/extract_proteins_from_augustus_output.py -i augustus.gff -o augustus.pep

~/soft/MAVR/scripts/hmmer3/parallel_hmmscan.py -i ~/data/Pfam/Pfam-A.hmm -s augustus.pep \
                                               -o augustus.pep.pfam.hits -c -t 10 \
                                               --domtblout augustus.pep.pfam.domtblout \
                                               --hmmer_dir ~/soft/hmmer-3.1b2-linux-intel-x86_64/binaries/

augustus_only_CDS.gff                                domtblout_dir       splited_fasta
augustus_only_CDS.gff.intersecting_with_repeats.gff  exe_list.t          tblout_dir
augustus_only_CDS.gff.intersecting_with_repeats.ids