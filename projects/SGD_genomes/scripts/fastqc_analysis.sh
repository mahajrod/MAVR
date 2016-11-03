#!/usr/bin/env bash

SAMPLE_LIST=($@)

for SAMPLE in ${SAMPLE_LIST[@]};
    do

    SAMPLE_GROUP=`echo ${SAMPLE} | cut -c1-4`

    mkdir -p ${FASTQC_STAT_DIR}/${SAMPLE_GROUP} ${FASTQC_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE};

    NUMBER_OF_FILES=`ls ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | wc -l`

    OUTPUT_DIR=${FASTQC_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE}/

    echo "Starting FastQC analysis"
    echo "    ${NUMBER_OF_FILES} files"

    FASTQC_STRING="${FASTQC_DIR}/fastqc -k 10 --nogroup -t ${NUMBER_OF_FILES} -o ${OUTPUT_DIR} ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/*"
    echo ${FASTQC_STRING}

    ${FASTQC_STRING}
    done