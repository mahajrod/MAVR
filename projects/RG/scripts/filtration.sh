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

    mkdir -p ${FILTERED_READS_DIR}/${SAMPLE_GROUP} ${FILTERED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE};
    mkdir -p ${FILTERING_STAT_DIR}/${SAMPLE_GROUP} ${FILTERING_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE};

    NUMBER_OF_FILES=`ls ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | wc -l`
    FILES=($(ls ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | sed 's/.gz//'));

    OUTPUT_STAT=${FILTERING_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE}/${SAMPLE}.filtering.stat
    OUTPUT_PREFIX=${FILTERED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/${SAMPLE}
    echo "Filtering reads by mean quality"
    echo "    ${NUMBER_OF_FILES} files"

    FACUT_STRING="${FACUT_BIN_DIR}/filter_by_mean_quality -t ${QUALITY_THRESHOLD} -q ${PHRED_SCORE_TYPE} -n ${READ_NAME_TYPE}  -f ${FILES[0]} -r ${FILES[1]} -p ${OUTPUT_PREFIX}"
    echo "${FACUT_STRING}  > ${OUTPUT_STAT}"

    ${FACUT_STRING}  > ${OUTPUT_STAT}
    done