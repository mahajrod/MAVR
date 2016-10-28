#!/usr/bin/env bash
#-----------------Environment variables-----------------
FASTQ_DIR=/home/genomerussia/main/fastq/
RAW_READS_DIR=${FASTQ_DIR}/raw/
UNPACKED_READS_DIR=${FASTQ_DIR}/unpacked/
FILTERED_READS_DIR=${FASTQ_DIR}/filtered/

ANALYSIS_DIR=/home/genomerussia/main/analysis/

ALIGNMENT_DIR=/home/genomerussia/main/analysis/alignment/
ALIGNMENT_BAM_DIR=${ALIGNMENT_DIR}/bam/
ALIGNMENT_LOG_DIR=${ALIGNMENT_DIR}/log/
ALIGNMENT_TMP_DIR=${ALIGNMENT_DIR}/tmp/

JF_DB_DIR=${ANALYSIS_DIR}/jf/

STAT_DIR=${ANALYSIS_DIR}/stat/
ADAPTERS_STAT_DIR=${STAT_DIR}/adapters/
FASTQC_STAT_DIR=${STAT_DIR}/fastqc/
FILTERING_STAT_DIR=${STAT_DIR}/filtering/
KMER_STAT_DIR=${STAT_DIR}/kmer/

TOOLS_DIR=/home/genomerussia/tools/
MAVR_SCRIPTS_DIR=/home/genomerussia/tools/MAVR/scripts/
FACUT_BIN_DIR=/home/genomerussia/tools/Facut/bin/
FASTQC_DIR=/home/genomerussia/tools/FastQC/
COOCKIECUTTER_SRC_DIR=/home/genomerussia/tools/Cookiecutter/src/

ADAPTER_KMER_FILE=${TOOLS_DIR}/service_sequences/trueseq_adapters_with_rev_com_23_mer.kmer
#PYTHONPATH=${PYTHONPATH}:/home/genomerussia/tools/MAVR
#export PYTONPATH
#-------------------------------------------------------

#----------------------Settings-------------------------
THREAD_NUMBER=60
KMER_SIZE=23
MEMORY=30G
PHRED_SCORE_TYPE=phred33
READ_NAME_TYPE=illumina
QUALITY_THRESHOLD=20
#-------------------------------------------------------


SAMPLE_LIST=($@)

for SAMPLE in ${SAMPLE_LIST[@]};
    do

    SAMPLE_GROUP=`echo ${SAMPLE} | cut -c1-4`

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

    LEFT_FILES_EXTRACT_STRING="zcat ${LEFT_FILES} > ${LEFT_UNPACKED_FILE} &"
    RIGHT_FILES_EXTRACT_STRING="zcat ${RIGHT_FILES} > ${RIGHT_UNPACKED_FILE}"

    echo "Unpacking and merging files for ${SAMPLE}"
    echo "    ${NUMBER_OF_FILES} files"

    echo ${LEFT_FILES_EXTRACT_STRING}
    ${LEFT_FILES_EXTRACT_STRING}

    echo ${RIGHT_FILES_EXTRACT_STRING}
    ${RIGHT_FILES_EXTRACT_STRING}

    done