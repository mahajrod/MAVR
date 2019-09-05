#!/usr/bin/env bash
#SBATCH --array=1-268%268
#SBATCH --time=99:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=16384       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=sol_homo
#SBbATCH --error=/work/toleksyk/toleksyk/errors/sol_homo.%A_%a.err
#SBATCH --output=/work/toleksyk/toleksyk/logs/sol_homo.%A_%a.out


PROFILE=/home/toleksyk/toleksyk/.profile
WORKDIR=/work/toleksyk/toleksyk/
SOFTDIR=${WORKDIR}/soft/
GENOMEDIR=${WORK}/genome/

EXONERATEDIR=${WORKDIR}/exonerate/
source ${PROFILE}

REFERENCE_SPECIES=homo_sapiens

SPLITED_PEP_DIR=${WORKDIR}/data/db/proteins/${REFERENCE_SPECIES}/ncbi/no_aligned/
PEP_FILE="${SPLITED_PEP_DIR}/no_aligned.${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.fasta"

SPECIES=solenodon_paradoxus
for AS in dovetail;
    do
    GENOME=${GENOMEDIR}/${SPECIES}/${AS}/${AS}.fasta
    EXONERATE_OUTPUT_DIR=${EXONERATEDIR}/${SPECIES}/${AS}/${REFERENCE_SPECIES}.noaligned/

    mkdir -p ${EXONERATE_OUTPUT_DIR}

    EXONERATE_OUTPUT_FILE=${EXONERATE_OUTPUT_DIR}/no_aligned.${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.exonerate.out


    exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}

    done


