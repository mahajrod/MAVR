#!/usr/bin/env bash
#SBATCH --array=1-948%948
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4096       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Bhol_eut
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/boechera_holboelli/logs/exonerate_eutrema.ncbi.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/boechera_holboelli/logs/exonerate_eutrema.ncbi.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

REFERENCE_SPECIES="eutrema_salsuginem"
SOURCE="ncbi"
SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/${REFERENCE_SPECIES}/${SOURCE}/pep/splited/"
PEP_FILE="${SPLITED_PEP_DIR}${REFERENCE_SPECIES}.${SOURCE}_${SLURM_ARRAY_TASK_ID}.fasta"

BHOLBOELLI_GENOMES_DIR="/work/pavlov/okochenova/genomes/boechera_holboelli/"
WORK_DIR="/work/pavlov/okochenova/genomes_workdir/boechera_holboelli/"

for SPECIES in SOAP_SOAP SSPACE_Sealer SSPACE_SOAP;
    do
    TARGET_GENOME="${BHOLBOELLI_GENOMES_DIR}${SPECIES}.fasta"
    TARGET_GENOME_DIR="${WORK_DIR}${SPECIES}/"
    EXONERATE_OUTPUT_DIR="${TARGET_GENOME_DIR}${REFERENCE_SPECIES}/"
    EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/exonerate_${REFERENCE_SPECIES}_pep_${SLURM_ARRAY_TASK_ID}.out"
    mkdir -p ${TARGET_GENOME_DIR}
    mkdir -p ${EXONERATE_OUTPUT_DIR}
    exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${TARGET_GENOME} > ${EXONERATE_OUTPUT_FILE}
    done


