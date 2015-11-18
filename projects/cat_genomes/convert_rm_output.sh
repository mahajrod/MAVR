#!/usr/bin/env bash

#supermicro

cd ~/data/genomes
for GENOME in *;
    do
    cd ${GENOME}/masking;
    #~/soft/MAVR/scripts/repeat_masking/convert_rm_out_to_gff.py -i final.assembly.fasta.out \
    #                                                            -o final.assembly.fasta.gff;
    sed -r 's/.*Class=(.*);Family.*/\1/' final.assembly.fasta.gff | sort | uniq > ${GENOME}_annotated_repeat_types.t
    cd ../../;
    done

