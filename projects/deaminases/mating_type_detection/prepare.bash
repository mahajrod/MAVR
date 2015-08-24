#!/usr/bin/env bash

source ./DESCRIPTION

cd ${WORKDIR}

mkdir -p ${KMER_DIR} ${MATING_RELATED_SEQUENCES_DIR} ${INDEX_DIR}

cd ${KMER_DIR}
~/Genetics/MAVR/scripts/kmer/get_kmer_list.py -i ${MATING_RELATED_SEQUENCES_DIR} -m 33 -a mating_rel_seqs -t 4 -r

BOWTIE2_DB_NAME="mating_rel_seqs"
bowtie2-build `ls -m ${MATING_RELATED_SEQUENCES_DIR}* | sed -r "s/, /,/g" | tr -d '\n'` ${INDEX_DIR}${BOWTIE2_DB_NAME}