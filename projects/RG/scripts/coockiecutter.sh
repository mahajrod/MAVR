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
#-------------------------------------------------------

SAMPLE_LIST=($@)

for SAMPLE in ${SAMPLE_LIST[@]};
    do

    SAMPLE_GROUP=`echo ${SAMPLE} | cut -c1-4`

    mkdir -p ${ADAPTERS_STAT_DIR}/${SAMPLE_GROUP} ${ADAPTERS_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE};

    NUMBER_OF_FILES=`ls ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | wc -l`
    FILES=($(ls ${UNPACKED_READS_DIR}/${SAMPLE_GROUP}/${SAMPLE}/* | sed 's/.gz//'));

    OUTPUT=${ADAPTERS_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE}/${SAMPLE}.adapters.stat

    echo "Counting reads with adapters"
    echo "    ${NUMBER_OF_FILES} files"

    COOCKIECUTTER_STRING="${COOCKIECUTTER_SRC_DIR}/counter -1 ${FILES[0]} -2 ${FILES[1]} -o ${ADAPTERS_STAT_DIR}/${SAMPLE_GROUP}/${SAMPLE}/ -f ${ADAPTER_KMER_FILE} > ${OUTPUT}"
    echo ${COOCKIECUTTER_STRING}

    ${COOCKIECUTTER_STRING}
    done