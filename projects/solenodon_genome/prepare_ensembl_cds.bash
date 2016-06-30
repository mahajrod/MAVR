#!/usr/bin/env bash
#supermicro

cd ~/data/Ensembl/release-70/fasta

for SPECIES in bos_taurus canis_familiaris equus_caballus homo_sapiens monodelphis_domestica mustela_putorius_furo mus_musculus;
    do

    mkdir -p  ${SPECIES}/cds_2;
    echo ${SPECIES};
    cd ${SPECIES}/cds_2;
    gffread ../../../gtf/${SPECIES}/${SPECIES}.gtf -g ../../${SPECIES}/genome/${SPECIES}.genome -x ${SPECIES}.cds -y ${SPECIES}.pep ;
    ~/soft/MAVR/scripts/annotation/get_transcript_to_pep_accordance_from_gtf.py -g ../../../gtf/${SPECIES}/${SPECIES}.gtf -o ${SPECIES}.transcript_to_pep_accordance;
    ~/soft/MAVR/scripts/annotation/trim_cds_and_remove_terminal_stop_codons.py -i ${SPECIES}.cds -o ${SPECIES}.trimmed.cds;
    ~/soft/MAVR/scripts/sequence/translate_sequences.py -i ${SPECIES}.trimmed.cds -o ${SPECIES}.trimmed.pep;
    ~/soft/MAVR/scripts/sequence/rename_sequence_ids.py -i ${SPECIES}.trimmed.pep -o ${SPECIES}.trimmed.renamed.pep -s ${SPECIES}.transcript_to_pep_accordance -k 0 -v 1 2>/dev/null;
    ~/soft/MAVR/scripts/sequence/compare_sequences.py -a ../pep/${SPECIES}.pep -b ${SPECIES}.trimmed.renamed.pep -o pep_comparison | tee ${SPECIES}.comparison.stat 2>/dev/null
    cd ../../;


    done