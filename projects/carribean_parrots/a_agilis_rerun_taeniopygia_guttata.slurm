#!/usr/bin/env bash
#SBATCH --array=1-80%80
#SBATCH --time=120:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=a_agi_taegut
#SBATCH --error=/work/toleksyk/toleksyk/errors/a_agi_taegut.%A_%a.err
#SBATCH --output=/work/toleksyk/toleksyk/logs/a_agi_taegut.%A_%a.out

PROFILE=/home/toleksyk/toleksyk/.profile
WORKDIR=/work/toleksyk/toleksyk/
SOFTDIR=${WORKDIR}/soft/
#EXONERATE_DIR
source ${PROFILE}

REFERENCE_SPECIES=taeniopygia_guttata

SPLITED_PEP_DIR=${WORKDIR}/data/db/proteins/${REFERENCE_SPECIES}/ncbi/tmp/agilis_rerun/splited/
PEP_FILE="${SPLITED_PEP_DIR}${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.fasta"


for SP in amazona_agilis;
    do
    GENOME=${WORK}/parrots/${SP}/genome_denovo/assemblies/10x/${SP}.10x.fasta
    EXONERATE_OUTPUT_DIR=${WORK}/parrots/${SP}/genome_denovo/annotation/exonerate/10x/${REFERENCE_SPECIES}/rerun/

    mkdir -p ${EXONERATE_OUTPUT_DIR}

    EXONERATE_OUTPUT_FILE=${EXONERATE_OUTPUT_DIR}${REFERENCE_SPECIES}_rerun_${SLURM_ARRAY_TASK_ID}.exonerate.out


    exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}

    done


