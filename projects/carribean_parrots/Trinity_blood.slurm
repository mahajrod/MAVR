#!/usr/bin/env bash
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=1192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Tri_blood
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --error=/work/toleksyk/toleksyk/errors/trinity_blood.%A_%a.err
#SBATCH --output=/work/toleksyk/toleksyk/logs/trinity_blood.%A_%a.out

WORKDIR=/work/toleksyk/toleksyk/
LOGDIR=${WORKDIR}/logs/
ERRORDIR=${WORKDIR}/errors/

DATADIR=${WORKDIR}/parrots/amazona_ventralis/transcriptome/reads/filtered/blood/trimmomatic/merged/

FORWARD=${DATADIR}/blood.pe_1.fastq
REVERSE=${DATADIR}/blood.pe_2.fastq


OUTDIR=${WORKDIR}/parrots/amazona_ventralis/transcriptome/assemblies/blood/
TRINITY_WORKDIR=${OUTPUT}/work/
mkdir -p ${OUTDIR}

SOFTDIR=${WORKDIR}/soft/

THREADS=8
GRID_NODE_THREADS=1
GRID_NODE_MAX_MEMORY=4G
GRID_EXEC=sbatch
GRID_CONF_FILE=


Trinity --seqType fq \
        --SS_lib_type RF \
        --left ${FORWARD} \
        --right ${REVERSE} \
        --output ${OUTPUT} \
        --workdir ${TRINITY_WORKDIR} \
        --CPU ${THREADS} \
        --grid_exec ${GRID_EXEC} \
        --grid_node_CPU ${GRID_NODE_THREADS} \
        --grid_node_max_memory ${GRID_NODE_MAX_MEMORY}
#--max_memory 20G