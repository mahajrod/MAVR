#!/usr/bin/env bash

WORKDIR=$1
TOOLS_DIR=/home/genomerussia/tools/

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

MAVR_SCRIPTS_DIR=${TOOLS_DIR}/MAVR/scripts/
FACUT_BIN_DIR=${TOOLS_DIR}/Facut/bin/
FASTQC_DIR=${TOOLS_DIR}/FastQC/
COOCKIECUTTER_SRC_DIR=${TOOLS_DIR}/Cookiecutter/src/

ADAPTER_KMER_FILE=${TOOLS_DIR}/service_sequences/trueseq_adapters_with_rev_com_23_mer.kmer

JF_DB_FILTERED_DIR=${ANALYSIS_DIR}/jf_filtered/
KMER_STAT_FILTERED_DIR=${STAT_DIR}/kmer_filtered/
FASTQC_STAT_FILTERED_DIR=${STAT_DIR}/fastqc_filtered/
FASTQC_STAT_RAW_SPLITED_DIR=${STAT_DIR}/fastqc_raw_splited/

#----------------------Settings-------------------------
THREAD_NUMBER=60
KMER_SIZE=23
MEMORY=30G
PHRED_SCORE_TYPE=phred33
READ_NAME_TYPE=illumina
QUALITY_THRESHOLD=20
#-------------------------------------------------------

echo "
WORKDIR=${WORKDIR}
TOOLS_DIR=/home/genomerussia/tools/

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

MAVR_SCRIPTS_DIR=${TOOLS_DIR}/MAVR/scripts/
FACUT_BIN_DIR=${TOOLS_DIR}/Facut/bin/
FASTQC_DIR=${TOOLS_DIR}/FastQC/
COOCKIECUTTER_SRC_DIR=${TOOLS_DIR}/Cookiecutter/src/

ADAPTER_KMER_FILE=${TOOLS_DIR}/service_sequences/trueseq_adapters_with_rev_com_23_mer.kmer

JF_DB_FILTERED_DIR=${ANALYSIS_DIR}/jf_filtered/
KMER_STAT_FILTERED_DIR=${STAT_DIR}/kmer_filtered/
FASTQC_STAT_FILTERED_DIR=${STAT_DIR}/fastqc_filtered/
FASTQC_STAT_RAW_SPLITED_DIR=${STAT_DIR}/fastqc_raw_splited/

#----------------------Settings-------------------------
THREAD_NUMBER=60
KMER_SIZE=23
MEMORY=30G
PHRED_SCORE_TYPE=phred33
READ_NAME_TYPE=illumina
QUALITY_THRESHOLD=20
#-------------------------------------------------------

export WORKDIR
export TOOLS_DIR
export FASTQ_DIR
export ANALYSIS_DIR
export ALIGNMENT_DIR
export RAW_READS_DIR
export UNPACKED_READS_DIR
export FILTERED_READS_DIR
export ALIGNMENT_BAM_DIR
export ALIGNMENT_LOG_DIR
export ALIGNMENT_TMP_DIR
export JF_DB_DIR
export STAT_DIR
export ADAPTERS_STAT_DIR
export FASTQC_STAT_DIR
export FILTERING_STAT_DIR
export KMER_STAT_DIR
export MAVR_SCRIPTS_DIR
export FACUT_BIN_DIR
export FASTQC_DIR
export COOCKIECUTTER_SRC_DIR
export ADAPTER_KMER_FILE

export THREAD_NUMBER
export KMER_SIZE
export MEMORY
export PHRED_SCORE_TYPE
export READ_NAME_TYPE
export QUALITY_THRESHOLD
#-------------------------------------------------------
" > $2

source $2
