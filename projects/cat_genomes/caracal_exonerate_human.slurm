#!/usr/bin/env bash
#SBATCH --array=1-890%333
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=cara_human
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/caracal/logs/exonerate_human.ncbi.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/caracal/logs/exonerate_human.ncbi.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/human/pep/human_ncbi_no_predicted/"
PEP_FILE="${SPLITED_PEP_DIR}human.ncbi.no_predicted_${SLURM_ARRAY_TASK_ID}.fasta"


CARACAL_GENOME="/work/pavlov/okochenova/genomes/caracal/final.assembly.fasta"

EXONERATE_OUTPUT_DIR="/work/pavlov/okochenova/genomes_workdir/caracal/exonerate_human_pep/"
mkdir -p ${EXONERATE_OUTPUT_DIR}
EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/exonerate_human_pep_${SLURM_ARRAY_TASK_ID}.out"


exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${CARACAL_GENOME} > ${EXONERATE_OUTPUT_FILE}