#!/usr/bin/env bash

THREADS=25
MAVR_DIR=~/soft/MAVR/
SCRIPTS_DIR=${MAVR_DIR}scripts/hmmer3/
#EXTRACT_HITS_SCRIPT=${SCRIPTS_DIR}extract_hits_from_hmm_output.py
HMMSCAN_SCRIPT=${SCRIPTS_DIR}parallel_hmmscan.py
TOP_HITS_SCRIPT=${SCRIPTS_DIR}extract_top_hits_from_hmm_output.py
FAMILIES_SCRIPT=${SCRIPTS_DIR}convert_top_hits_to_families.py

TREEFAM_HMM3_PROFILE=~/data/TreeFam_new/treefam9.hmm3

for SPECIES in pantholops_hodgsonii hippotragus_niger ovis_aries;
    do
    cd ~/workdir/sable/${SPECIES}/
    $HMMSCAN_SCRIPT -t ${THREADS} \
                    -s ${SPECIES}.pep \
                    -i ${TREEFAM_HMM3_PROFILE} \
                    -c -o ${SPECIES}_hmmscan.report \
                    -d hmmscan_output_dir/

    $TOP_HITS_SCRIPT -t ${THREADS} \
                     -r -d top_hits_dir/ \
                     -i hmmscan_output_dir/ \
                     -o ${SPECIES}_hmmscan_top_hits.tab \
                     -f hmmer3-text -n ${SPECIES}_hmmscan_not_found.t \
                     -g ${SPECIES}_hmmscan_not_significant.t

    $FAMILIES_SCRIPT -i ${SPECIES}_hmmscan_top_hits.tab \
                    -e -k 1 -a 0 \
                    -o ${SPECIES}_hmmscan_families.tab
    done

