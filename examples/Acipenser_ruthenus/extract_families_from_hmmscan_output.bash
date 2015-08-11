#!/usr/bin/env bash

#server=dell

MAVR_DIR=~/Soft/MAVR/
SCRIPTS_DIR=${MAVR_DIR}scripts/hmmer3/
#EXTRACT_HITS_SCRIPT=${SCRIPTS_DIR}extract_hits_from_hmm_output.py
TOP_HITS_SCRIPT=${SCRIPTS_DIR}extract_top_hits_from_hmm_output.py
FAMILIES_SCRIPT=${SCRIPTS_DIR}convert_top_hits_to_families.py
SPECIES="acipenser_ruthenus"
cd ~/workdir/Acipenser_ruthenus/transcriptomes/RNA1/TransDecoder

$TOP_HITS_SCRIPT -i ~/workdir/Acipenser_ruthenus/transcriptomes/RNA1/TransDecoder/splited_output_dir/ \
                    -o ${SPECIES}_hmmscan_top_hits.tab \
                    -f hmmer3-text -n ${SPECIES}_hmmscan_not_found.t \
                    -g ${SPECIES}_hmmscan_not_significant.t

$FAMILIES_SCRIPT -i ${SPECIES}_hmmscan_top_hits.tab \
                    -e -k 1 -a 0 \
                    -o ${SPECIES}_hmmscan_families.tab



