#!/usr/bin/env bash
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Avent_trans_blood
#SBATCH --ntasks=8
#SBATCH --error=/work/toleksyk/toleksyk/errors/trimmomatic_blood.%A_%a.err
#SBATCH --output/work/toleksyk/toleksyk/logs/trimmomatic_blood.%A_%a.out

WORKDIR=/work/toleksyk/toleksyk/
LOGDIR=${WORKDIR}/logs/
ERRORDIR=${WORKDIR}/errors/

DATADIR=${WORKDIR}/parrots/amazona_ventralis/transcriptome/reads/filtered/blood/stirka/
OUTDIR=${WORKDIR}/parrots/amazona_ventralis/transcriptome/reads/filtered/blood/trimmomatic/

SOFTDIR=${WORKDIR}/soft/

module load trimmomatic/0.36

THREADS=8
TRIMMOMATIC_ADAPTERS=${SOFT}/Trimmomatic-0.38/adapters/TruSeq2-PE.fa

mkdir -p ${OUTDIR}

for SAMPLE in Parrot140  Parrot335  Parrot341;
    do
    mkdir -p ${OUTDIR}/${SAMPLE}/

    FORWARD=${DATADIR}/${SAMPLE}/${SAMPLE}.trim_1.fastq.gz
    REVERSE=${DATADIR}/${SAMPLE}/${SAMPLE}.trim_2.fastq.gz

    OUTPUT_PREFIX=${OUTDIR}/${SAMPLE}/${SAMPLE}

    FORWARD_PE_OUT=${OUTPUT_PREFIX}.pe_1.fastq
    FORWARD_SE_OUT=${OUTPUT_PREFIX}.se_1.fastq
    FORWARD_PE_OUT=${OUTPUT_PREFIX}.pe_2.fastq
    FORWARD_SE_OUT=${OUTPUT_PREFIX}.se_2.fastq

    TRIMMOMATIC_LOG=${OUTPUT_PREFIX}.trimmomatic.log

    trimmomatic PE -threads ${THREADS} -phred33 ${FORWARD} ${REVERSE} ${FORWARD_PE_OUT} ${FORWARD_SE_OUT} ${REVERSE_PE_OUT} ${REVERSE_SE_OUT} ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:2:30:10:1 SLIDINGWINDOW:20:20 MINLEN:50 > ${TRIMMOMATIC_LOG} 2>&1

    done

for SAMPLE in Parrot13;
    do
    mkdir -p ${OUTDIR}/${SAMPLE}/

    FORWARD=${DATADIR}/${SAMPLE}/${SAMPLE}.trim_1.fastq
    REVERSE=${DATADIR}/${SAMPLE}/${SAMPLE}.trim_2.fastq

    OUTPUT_PREFIX=${OUTDIR}/${SAMPLE}/${SAMPLE}

    FORWARD_PE_OUT=${OUTPUT_PREFIX}.pe_1.fastq
    FORWARD_SE_OUT=${OUTPUT_PREFIX}.se_1.fastq
    FORWARD_PE_OUT=${OUTPUT_PREFIX}.pe_2.fastq
    FORWARD_SE_OUT=${OUTPUT_PREFIX}.se_2.fastq

    TRIMMOMATIC_LOG=${OUTPUT_PREFIX}.trimmomatic.log

    trimmomatic PE -threads ${THREADS} -phred33 ${FORWARD} ${REVERSE} ${FORWARD_PE_OUT} ${FORWARD_SE_OUT} ${REVERSE_PE_OUT} ${REVERSE_SE_OUT} ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:2:30:10:1 SLIDINGWINDOW:20:20 MINLEN:50 > ${TRIMMOMATIC_LOG} 2>&1

    done

DATADIR2=${WORKDIR}/parrots/amazona_ventralis/transcriptome/reads/filtered/liver/stirka/
OUTDIR2=${WORKDIR}/parrots/amazona_ventralis/transcriptome/reads/filtered/liver/trimmomatic/

for SAMPLE in Parrot13;
    do
    mkdir -p ${OUTDIR2}/${SAMPLE}/

    FORWARD=${DATADIR2}/${SAMPLE}/${SAMPLE}.trim_1.fastq.gz
    REVERSE=${DATADIR2}/${SAMPLE}/${SAMPLE}.trim_2.fastq.gz

    OUTPUT_PREFIX=${OUTDIR2}/${SAMPLE}/${SAMPLE}

    FORWARD_PE_OUT=${OUTPUT_PREFIX}.pe_1.fastq
    FORWARD_SE_OUT=${OUTPUT_PREFIX}.se_1.fastq
    FORWARD_PE_OUT=${OUTPUT_PREFIX}.pe_2.fastq
    FORWARD_SE_OUT=${OUTPUT_PREFIX}.se_2.fastq

    TRIMMOMATIC_LOG=${OUTPUT_PREFIX}.trimmomatic.log

    trimmomatic PE -threads ${THREADS} -phred33 ${FORWARD} ${REVERSE} ${FORWARD_PE_OUT} ${FORWARD_SE_OUT} ${REVERSE_PE_OUT} ${REVERSE_SE_OUT} ILLUMINACLIP:${TRIMMOMATIC_ADAPTERS}:2:30:10:1 SLIDINGWINDOW:20:20 MINLEN:50 > ${TRIMMOMATIC_LOG} 2>&1

    done