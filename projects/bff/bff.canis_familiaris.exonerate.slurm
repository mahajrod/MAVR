#!/usr/bin/env bash
#SBATCH --array=1-980%500
#SBATCH --time=120:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=6192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=bff_can
#SBATCH --error=/work/toleksyk/toleksyk/tmp/mustela_nigripes/log/bff_can.%A_%a.err
#SBATCH --output=/work/toleksyk/toleksyk/tmp/mustela_nigripes/log/bff_can.%A_%a.out



PROFILE=/home/toleksyk/toleksyk/.profile
WORKDIR=/work/toleksyk/toleksyk/
SOFTDIR=${WORKDIR}/soft/
#EXONERATE_DIR
source ${PROFILE}

REFERENCE_SPECIES=canis_familiaris

SPLITED_PEP_DIR=${WORKDIR}/data/db/proteins/${REFERENCE_SPECIES}/ncbi/splited_pep/
PEP_FILE="${SPLITED_PEP_DIR}${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.fasta"


for SP in amazona_vittata;
    do
    GENOME=${WORK}/tmp/mustela_nigripes/mustela_nigripes.v2.smithsonian.fasta
    EXONERATE_OUTPUT_DIR=${WORK}/tmp/mustela_nigripes/exonerate/${REFERENCE_SPECIES}/

    mkdir -p ${EXONERATE_OUTPUT_DIR}

    EXONERATE_OUTPUT_FILE=${EXONERATE_OUTPUT_DIR}${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.exonerate.out


    exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}

    done


