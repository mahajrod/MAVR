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