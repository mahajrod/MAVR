#!/usr/bin/env bash
#SBATCH --array=723
#SBATCH --time=120:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=16192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=sol_con
#SBATCH --error=/work/toleksyk/toleksyk/errors/sol_con.%A_%a.err
#SBATCH --output=/work/toleksyk/toleksyk/logs/sol_con.%A_%a.out


PROFILE=/home/toleksyk/toleksyk/.profile
WORKDIR=/work/toleksyk/toleksyk/
SOFTDIR=${WORKDIR}/soft/
GENOMEDIR=${WORK}/genome/

EXONERATEDIR=${WORKDIR}/exonerate/
source ${PROFILE}

REFERENCE_SPECIES=condylura_cristata

SPLITED_PEP_DIR=${WORKDIR}/data/db/proteins/${REFERENCE_SPECIES}/ncbi/ConCri1.0/splited_pep/
PEP_FILE="${SPLITED_PEP_DIR}${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.fasta"

SPECIES=solenodon_paradoxus
for AS in dovetail;
    do
    GENOME=${GENOMEDIR}/${SPECIES}/${AS}/${AS}.fasta
    EXONERATE_OUTPUT_DIR=${EXONERATEDIR}/${SPECIES}/${AS}/${REFERENCE_SPECIES}/

    mkdir -p ${EXONERATE_OUTPUT_DIR}

    EXONERATE_OUTPUT_FILE=${EXONERATE_OUTPUT_DIR}/${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.exonerate.rerun2.out


    exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}

    done


