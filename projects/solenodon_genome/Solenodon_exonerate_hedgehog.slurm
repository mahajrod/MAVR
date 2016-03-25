#!/usr/bin/env bash
#SBATCH --array=1-935%935
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=10000       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Sol_hedge
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/solenodon/logs/exonerate_hedgehog.ncbi.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/solenodon/logs/exonerate_hedgehog.ncbi.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/erinaceus_europaeus/ncbi/pep/splited/"
PEP_FILE="${SPLITED_PEP_DIR}erinaceus_europaeus.ncbi_${SLURM_ARRAY_TASK_ID}.fasta"

GENOME="/work/pavlov/okochenova/genomes/solenodon/scaffolds.fa"

EXONERATE_OUTPUT_DIR="/work/pavlov/okochenova/genomes_workdir/solenodon/exonerate_hedgehog_pep/"
mkdir -p ${EXONERATE_OUTPUT_DIR}
EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/exonerate_hedgehog_pep_${SLURM_ARRAY_TASK_ID}.out"

exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}
