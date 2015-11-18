#!/usr/bin/env bash

cd ~/data/genomes
for GENOME in *;
    do
    cd ${GENOME}/masking;

    grep "Class=DNA|Class=DNA?|Class=LINE|Class=LTR|Class=LTR?|Class=RC|Class=RC?|Class=SINE|Class=SINE?" final.assembly.fasta.gff > final.assembly.fasta_selected_classes.gff

    cd ../../;
    done