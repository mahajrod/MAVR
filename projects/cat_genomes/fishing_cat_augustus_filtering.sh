#!/usr/bin/env bash

cd /home/skliver/workdir/cats/fishing_cat

grep -P "\t\CDSt" fishing_cat_augustus_cat_model.gff > fishing_cat_augustus_cat_model_CDS_only.gff
bedtools intersect -u -a fishing_cat_augustus_cat_model_CDS_only.gff \
                   -b ~/data/genomes/fishing_cat/masking/final.assembly.fasta_selected_classes.gff > fishing_cat_CDS_intersecting_with_repeats.gff


sed "s/.*=//" fishing_cat_CDS_intersecting_with_repeats.gff | sort | uniq > fishing_cat_genes_with_CDS_intersecting_with_repeats.ids