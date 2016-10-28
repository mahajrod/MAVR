#!/usr/bin/env bash

#for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123; do  mkdir -p ../../unpacked/GR00/${SAMPLE}; array=($(ls ${SAMPLE}/ | sed 's/.gz//')); echo ${array[0]};  gunzip -c ${SAMPLE}/${array[0]}.gz > ../../unpacked/GR00/${SAMPLE}/${array[0]} &  echo ${array[1]}; gunzip -c ${SAMPLE}/${array[1]}.gz > ../../unpacked/GR00/${SAMPLE}/${array[1]};  done

for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123;
    do
    mkdir -p ../../unpacked/GR00/${SAMPLE};
    FILES=($(ls ${SAMPLE}/ | sed 's/.gz//'));
    echo ${FILES[0]};
    gunzip -c ${SAMPLE}/${FILES[0]}.gz > ../../unpacked/GR00/${SAMPLE}/${FILES[0]} &
    echo ${FILES[1]};
    gunzip -c ${SAMPLE}/${FILES[1]}.gz > ../../unpacked/GR00/${SAMPLE}/${FILES[1]};
    done

RAW_READS_DIR=/home/genomerussia/main/fastq/raw/
UNPACKED_READS_DIR=/home/genomerussia/main/fastq/unpacked/

for SAMPLE in GR0331 GR0332 GR0333 GR0334 GR0335 GR0336 GR0338 GR0339 GR0340 GR0383 GR0384;
    do
    SAMPLE_GROUP=`echo ${SAMPLE} | cut -c1-4`

    mkdir -p ${UNPACKED_READS_DIR}/${SAMPLE_GROUP} ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE};
    FILES=($(ls ${SAMPLE}/ | sed 's/.gz//'));
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