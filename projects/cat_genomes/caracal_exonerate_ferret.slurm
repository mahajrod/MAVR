#!/usr/bin/env bash
#SBATCH --array=1-995%333
#SBATCH --time=100:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=8192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=cara_ferret
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/caracal/logs/exonerate_ferret.ncbi.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/caracal/logs/exonerate_ferret.ncbi.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/ferret/pep/splited_ncbi_pep/"

PEP_FILE="${SPLITED_PEP_DIR}mustela_putorius_furo_${SLURM_ARRAY_TASK_ID}.fasta"
CARACAL_GENOME="/work/pavlov/okochenova/genomes/caracal/final.assembly.fasta"

mkdir -p /work/pavlov/okochenova/genomes_workdir/caracal/exonerate_ferret_ncbi_pep/
EXONERATE_OUTPUT_FILE="/work/pavlov/okochenova/genomes_workdir/caracal/exonerate_ferret_ncbi_pep/exonerate_ferret_ncbi_pep_${SLURM_ARRAY_TASK_ID}.out"

exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${CARACAL_GENOME} > ${EXONERATE_OUTPUT_FILE}