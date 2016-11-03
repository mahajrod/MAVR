#!/usr/bin/env bash


SAMPLE_LIST=($@)

for SAMPLE in ${SAMPLE_LIST[@]};
    do

    SAMPLE_GROUP=`echo ${SAMPLE} | cut -c1-9`

    mkdir -p ${ADAPTERS_STAT_DIR}/${SAMPLE_GROUP} ${ADAPTERS_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE};

    NUMBER_OF_FILES=`ls ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | wc -l`
    FILES=($(ls ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | sed 's/.gz//'));

    OUTPUT=${ADAPTERS_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE}/${SAMPLE}.adapters.stat

    echo "Counting reads with adapters"
    echo "    ${NUMBER_OF_FILES} files"

    COOCKIECUTTER_STRING="${COOCKIECUTTER_SRC_DIR}/counter -1 ${FILES[0]} -2 ${FILES[1]} -o ${ADAPTERS_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE}/ -f ${ADAPTER_KMER_FILE}"
    echo "${COOCKIECUTTER_STRING} > ${OUTPUT}"

    ${COOCKIECUTTER_STRING} > ${OUTPUT}
    done