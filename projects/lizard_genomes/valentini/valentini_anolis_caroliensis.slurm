#!/usr/bin/env bash
#SBATCH --array=1-968%400
#SBATCH --time=120:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=8192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=val_anol
#SBATCH --error=/work/toleksyk/toleksyk/errors/val_anol.%A_%a.err
#SBATCH --output=/work/toleksyk/toleksyk/logs/val_anol.%A_%a.out


PROFILE=/home/toleksyk/toleksyk/.profile
WORKDIR=/work/toleksyk/toleksyk/
SOFTDIR=${WORKDIR}/soft/
GENOMEDIR=${WORK}/genome/

EXONERATEDIR=${WORKDIR}/exonerate/
source ${PROFILE}

REFERENCE_SPECIES=anolis_carolinensis

SPLITED_PEP_DIR=${WORKDIR}/data/db/proteins/${REFERENCE_SPECIES}/ncbi/splited/
PEP_FILE="${SPLITED_PEP_DIR}${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.fasta"


for SP in darevskia_valentini;
    do
    GENOME=${GENOMEDIR}/darevskia_valentini/darevskia_valentini.v2.fasta
    EXONERATE_OUTPUT_DIR=${EXONERATEDIR}/${SP}/${REFERENCE_SPECIES}/

    mkdir -p ${EXONERATE_OUTPUT_DIR}

    EXONERATE_OUTPUT_FILE=${EXONERATE_OUTPUT_DIR}/${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.exonerate.out


    exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}

    done


