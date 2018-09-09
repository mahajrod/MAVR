#!/usr/bin/env bash
#SBATCH --array=1-973%400
#SBATCH --time=120:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=A_vit_t_gut
#SBATCH --error=/work/toleksyk/toleksyk/errors/A_vit_t_gut.%A_%a.err
#SBATCH --output=/work/toleksyk/toleksyk/logs/A_vit_t_gut.%A_%a.out

PROFILE=/home/toleksyk/toleksyk/.profile
WORKDIR=/work/toleksyk/toleksyk/
SOFTDIR=${WORKDIR}/soft/
#EXONERATE_DIR
source ${PROFILE}

REFERENCE_SPECIES=taeniopygia_guttata

SPLITED_PEP_DIR=${WORKDIR}/data/db/proteins/${REFERENCE_SPECIES}/ncbi/splited/
PEP_FILE="${SPLITED_PEP_DIR}${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.fasta"


for SP in amazona_vittata;
    do
    GENOME=${WORK}/parrots/${SP}/genome_denovo/assemblies/fermi_SSPACE_gapcloser_LRNA/L_RNA_scaffolder.fasta
    EXONERATE_OUTPUT_DIR=${WORK}/parrots/${SP}/genome_denovo/annotation/exonerate/fermi_SSPACE_gapcloser_LRNA/${REFERENCE_SPECIES}/

    mkdir -p ${EXONERATE_OUTPUT_DIR}

    EXONERATE_OUTPUT_FILE=${EXONERATE_OUTPUT_DIR}${REFERENCE_SPECIES}_${SLURM_ARRAY_TASK_ID}.exonerate.out


    exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}

    done


