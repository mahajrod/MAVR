#!/usr/bin/env bash

# supermicro

MAVR_DIR=~/soft/MAVR/
SCRIPTS_DIR=${MAVR_DIR}scripts/expansion_hmm/
#EXTRACT_HITS_SCRIPT=${SCRIPTS_DIR}extract_hits_from_hmm_output.py
TOP_HITS_SCRIPT=${SCRIPTS_DIR}extract_top_hits_from_hmm_output.py
FAMILIES_SCRIPT=${SCRIPTS_DIR}convert_top_hits.py

for SPECIES in sable;
    do
    cd ~/workdir/sable/
    $TOP_HITS_SCRIPT -i ${SPECIES}_hmmscan.report \
                    -o ${SPECIES}_hmmscan_top_hits.tab \
                    -f hmmer3-text -n ${SPECIES}_hmmscan_not_found.t \
                    -g ${SPECIES}_hmmscan_not_significant.t

    $FAMILIES_SCRIPT -i ${SPECIES}_hmmscan_top_hits.tab \
                    -e -k 1 -a 0 \
                    -o ${SPECIES}_hmmscan_families.tab
    done


