#!/usr/bin/env bash

for SAMPLE in *; do cd ${SAMPLE}; ~/Soft/MAVR/scripts/annotation/parallel_augustus.py -i ${SAMPLE}.fixed.masked_selected_repeat_types.fa -o ${SAMPLE}.augustus -t 32 -x ${SAMPLE%_*} -s human -c ~/Soft/augustus-3.2.1/config/ -p /mnt/guatemala/skliver/data/Pfam/Pfam-A.hmm -w /mnt/guatemala/skliver/data/SwissProt/swissprot --softmasking --hintsfile ${SAMPLE}.all_hints.gff --extrinsicCfgFile ~/Soft/augustus-3.2.1/config/extrinsic/extrinsic.RM.EXNT.EXNS.cfg -a ~/Soft/augustus-3.2.1/bin/; cd ../; done
