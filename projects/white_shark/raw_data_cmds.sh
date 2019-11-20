#!/usr/bin/env bash


skliver@dell:/mnt/guatemala/skliver/white_shark_project/raw_data$ basename -s _alignment.fasta `ls monoclusters_pep_alignment` | xargs -P 15 -I NAME ~/Soft/MAVR/scripts/multiple_alignment/get_codon_alignment.py -p monoclusters_pep_alignment/NAME_alignment.fasta -c labeled_cds_fixed_ids/ -o monoclusters_cds_alignment/NAME_alignment.fasta -a merged.labeled.accordance -r -i nuc_tmp.idx


