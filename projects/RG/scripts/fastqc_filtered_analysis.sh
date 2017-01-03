#!/usr/bin/env bash
#-----------------Environment variables-----------------
WORKDIR=/home/genomerussia/main/
FASTQ_DIR=${WORKDIR}/fastq/
ANALYSIS_DIR=${WORKDIR}/analysis/
ALIGNMENT_DIR=${WORKDIR}/analysis/alignment/

RAW_READS_DIR=${FASTQ_DIR}/raw/
UNPACKED_READS_DIR=${FASTQ_DIR}/unpacked/
FILTERED_READS_DIR=${FASTQ_DIR}/filtered/

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
MAVR_SCRIPTS_DIR=${TOOLS_DIR}/MAVR/scripts/
FACUT_BIN_DIR=${TOOLS_DIR}/Facut/bin/
FASTQC_DIR=${TOOLS_DIR}/FastQC/
COOCKIECUTTER_SRC_DIR=${TOOLS_DIR}/Cookiecutter/src/

ADAPTER_KMER_FILE=${TOOLS_DIR}/service_sequences/trueseq_adapters_with_rev_com_23_mer.kmer
#PYTHONPATH=${PYTHONPATH}:/home/genomerussia/tools/MAVR
#export PYTONPATH

JF_DB_FILTERED_DIR=${ANALYSIS_DIR}/jf_filtered/
KMER_STAT_FILTERED_DIR=${STAT_DIR}/kmer_filtered/
FASTQC_STAT_FILTERED_DIR=${STAT_DIR}/fastqc_filtered/
FASTQC_STAT_RAW_SPLITED_DIR=${STAT_DIR}/fastqc_raw_splited/

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

    mkdir -p ${FASTQC_STAT_FILTERED_DIR}/${SAMPLE_GROUP} ${FASTQC_STAT_FILTERED_DIR}/${SAMPLE_GROUP}/${SAMPLE};

    NUMBER_OF_FILES=`ls ${FILTERED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | wc -l`

    OUTPUT_DIR=${FASTQC_STAT_FILTERED_DIR}/${SAMPLE_GROUP}/${SAMPLE}/

    echo "Starting FastQC analysis"
    echo "    ${NUMBER_OF_FILES} files"

    FASTQC_STRING="${FASTQC_DIR}/fastqc -k 10 --nogroup -t ${NUMBER_OF_FILES} -o ${OUTPUT_DIR} ${FILTERED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/*"
    echo ${FASTQC_STRING}

    ${FASTQC_STRING}
    done