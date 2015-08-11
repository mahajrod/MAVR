#!/usr/bin/env bash

# supermicro

MAVR_DIR=~/soft/MAVR/
SCRIPTS_DIR=${MAVR_DIR}scripts/hmmer3/
#EXTRACT_HITS_SCRIPT=${SCRIPTS_DIR}extract_hits_from_hmm_output.py
TOP_HITS_SCRIPT=${SCRIPTS_DIR}extract_top_hits_from_hmm_output.py
FAMILIES_SCRIPT=${SCRIPTS_DIR}convert_top_hits_to_families.py
SPECIES="mammalia"
cd ~/data/TreeFam_new/check_ambigious/

$TOP_HITS_SCRIPT -i ~/data/TreeFam_new/check_ambigious/splited_output_dir/ \
                    -o ${SPECIES}_hmmscan_top_hits.tab \
                    -f hmmer3-text -n ${SPECIES}_hmmscan_not_found.t \
                    -g ${SPECIES}_hmmscan_not_significant.t \
                    -t 20

$FAMILIES_SCRIPT -i ${SPECIES}_hmmscan_top_hits.tab \
                    -e -k 1 -a 0 \
                    -o ${SPECIES}_hmmscan_families.tab



