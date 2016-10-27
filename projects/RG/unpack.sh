#!/usr/bin/env bash

SAMPLE_LIST=($@)

RAW_READS_DIR=/home/genomerussia/main/fastq/raw/
UNPACKED_READS_DIR=/home/genomerussia/main/fastq/unpacked/

for SAMPLE in ${SAMPLE_LIST[@]};
    do
    echo "Unpacking and merging files for ${SAMPLE}"
    SAMPLE_GROUP=`echo ${SAMPLE} | cut -c1-4`

    mkdir -p ${UNPACKED_READS_DIR}/${SAMPLE_GROUP} ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE};
    FILES=($(ls ${RAW_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/ | sed 's/.gz//'));
    NUMBER_OF_FILES=${#FILES[@]}
    LEFT_FILES=
    RIGHT_FILES=
    INDEX=0

    echo ${INDEX}
    echo ${NUMBER_OF_FILES}
    while ((INDEX < NUMBER_OF_FILES));
        do
        LEFT_FILES="${LEFT_FILES} ${FILES[${INDEX}]}"
        RIGHT_FILES="${RIGHT_FILES} ${FILES[`expr ${INDEX}+1`]}"
        INDEX=`expr $INDEX + 2`

        done

    echo ${LEFT_FILES}
    echo ${RIGHT_FILES}

    zcat ${LEFT_FILES} > ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/${SAMPLE}_1.fastq &
    zcat ${RIGHT_FILES} > ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/${SAMPLE}_2.fastq

    done