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

#PYTHONPATH=${PYTHONPATH}:/home/genomerussia/tools/MAVR
#export PYTONPATH
#-------------------------------------------------------

#----------------------Settings-------------------------

THREAD_NUMBER=60
KMER_SIZE=23
MEMORY=30G
#-------------------------------------------------------

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