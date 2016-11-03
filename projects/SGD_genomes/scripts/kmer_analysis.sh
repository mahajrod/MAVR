#!/usr/bin/env bash


SAMPLE_LIST=($@)

for SAMPLE in ${SAMPLE_LIST[@]};
    do

    SAMPLE_GROUP=`echo ${SAMPLE} | cut -c1-4`
    SAMPLE_JF_DIR=${JF_DB_DIR}/${SAMPLE_GROUP}/${SAMPLE}/
    SAMPLE_KMER_STAT_DIR=${KMER_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE}/

    mkdir -p ${KMER_STAT_DIR}/${SAMPLE_GROUP} ${SAMPLE_KMER_STAT_DIR};
    mkdir -p ${JF_DB_DIR}/${SAMPLE_GROUP} ${SAMPLE_JF_DIR};
    #get comma-separated list of files in ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}
    FILES_COMMA=`ls -m ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | sed -r "s/, /,/g" | tr -d '\n'`;
    NUMBER_OF_FILES=`ls ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | wc -l`

    OUTPUT_PREFIX=${SAMPLE_JF_DIR}/${SAMPLE}

    echo "Counting k-mer distribution for ${SAMPLE}"
    echo "    ${NUMBER_OF_FILES} files"

    PYTHONPATH=${PYTHONPATH}:/home/genomerussia/tools/MAVR
    KMER_STRING=" ${MAVR_SCRIPTS_DIR}/kmer/draw_kmer_distribution_from_fastq.py -i ${FILES_COMMA} -t ${THREAD_NUMBER} -m ${KMER_SIZE} -b -s ${MEMORY} -e png -w 3 -g 80 -o ${OUTPUT_PREFIX}"
    echo ${KMER_STRING}

    ${KMER_STRING}

    echo "Coping statistics to statistics directory"
    CP_STRING="cp ${SAMPLE_JF_DIR}/${SAMPLE}_${KMER_SIZE}_mer.histo ${SAMPLE_JF_DIR}/${SAMPLE}_${KMER_SIZE}_mer_histo* ${SAMPLE_KMER_STAT_DIR}"
    echo ${CP_STRING}
    ${CP_STRING}

    done