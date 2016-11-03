#!/usr/bin/env bash

SAMPLE_LIST=($@)

for SAMPLE in ${SAMPLE_LIST[@]};
    do

    SAMPLE_GROUP=`echo ${SAMPLE} | cut -c1-9`

    mkdir -p ${UNPACKED_READS_DIR}/${SAMPLE_GROUP} ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE};
    FILES=($(ls ${RAW_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | sed 's/.gz//'));

    NUMBER_OF_FILES=${#FILES[@]}
    LEFT_FILES=
    RIGHT_FILES=
    INDEX=0

    #echo ${INDEX}
    #echo ${NUMBER_OF_FILES}
    while ((INDEX < NUMBER_OF_FILES));
        do
        LEFT_FILES="${LEFT_FILES} ${FILES[${INDEX}]}"
        RIGHT_FILES="${RIGHT_FILES} ${FILES[`expr ${INDEX}+1`]}"
        INDEX=`expr $INDEX + 2`

        done

    #echo ${LEFT_FILES}
    #echo ${RIGHT_FILES}

    LEFT_UNPACKED_FILE=${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/${SAMPLE}_1.fastq
    RIGHT_UNPACKED_FILE=${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/${SAMPLE}_2.fastq

    LEFT_FILES_EXTRACT_STRING="zcat ${LEFT_FILES}"
    RIGHT_FILES_EXTRACT_STRING="zcat ${RIGHT_FILES}"

    echo "Unpacking and merging files for ${SAMPLE}"
    echo "    ${NUMBER_OF_FILES} files"

    echo "${LEFT_FILES_EXTRACT_STRING} > ${LEFT_UNPACKED_FILE} &"
    ${LEFT_FILES_EXTRACT_STRING} > ${LEFT_UNPACKED_FILE} &

    echo "${RIGHT_FILES_EXTRACT_STRING} > ${RIGHT_UNPACKED_FILE}"
    ${RIGHT_FILES_EXTRACT_STRING} > ${RIGHT_UNPACKED_FILE}

    done