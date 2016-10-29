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

    #${KMER_STRING}

    cp ${SAMPLE_JF_DIR}/${SAMPLE}_${KMER_SIZE}_mer.histo ${SAMPLE_JF_DIR}/${SAMPLE}_${KMER_SIZE}_mer_histo* ${SAMPLE_KMER_STAT_DIR}
    done