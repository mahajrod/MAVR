#!/usr/bin/env bash

SAMPLE_LIST=($@)

for SAMPLE in ${SAMPLE_LIST[@]};
    do

    #SAMPLE_GROUP=`echo ${SAMPLE} | cut -c1-4`

    mkdir -p ${FILTERED_READS_DIR} ${FILTERED_READS_DIR}/${SAMPLE};
    mkdir -p ${FILTERING_STAT_DIR} ${FILTERING_STAT_DIR}/${SAMPLE};

    NUMBER_OF_FILES=`ls ${UNPACKED_READS_DIR}/${SAMPLE}/* | wc -l`
    FILES=($(ls ${UNPACKED_READS_DIR}/${SAMPLE}/* | sed 's/.gz//'));

    OUTPUT_STAT=${FILTERING_STAT_DIR}/${SAMPLE}/${SAMPLE}.filtering.stat
    OUTPUT_PREFIX=${FILTERED_READS_DIR}/${SAMPLE}/${SAMPLE}
    echo "Filtering reads by mean quality"
    echo "    ${NUMBER_OF_FILES} files"

    FACUT_STRING="${FACUT_BIN_DIR}/filter_by_mean_quality -t ${QUALITY_THRESHOLD} -q ${PHRED_SCORE_TYPE} -n ${READ_NAME_TYPE}  -f ${FILES[0]} -r ${FILES[1]} -p ${OUTPUT_PREFIX}"
    echo "${FACUT_STRING}  > ${OUTPUT_STAT}"

    ${FACUT_STRING}  > ${OUTPUT_STAT}
    done