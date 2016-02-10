#!/usr/bin/env bash
#SBATCH --array=1-900%900
#SBATCH --time=100:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=8192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=leocat_ferret
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/leopard_cat/logs/exonerate_ferret.ncbi.replaced_selenocystein_%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/leopard_cat/logs/exonerate_ferret.ncbi.replaced_selenocystein_%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/ferret/pep/splited_ncbi_pep_replaced_selenocystein/"

PEP_FILE="${SPLITED_PEP_DIR}mustela_putorius_furo.ncbi.replaced_selenocystein_${SLURM_ARRAY_TASK_ID}.fasta"
GENOME="/work/pavlov/okochenova/genomes/leopard_cat/final.assembly.fasta"

OUT_DIR="/work/pavlov/okochenova/genomes_workdir/leopard_cat/exonerate_ferret_ncbi_pep/"
mkdir -p ${OUT_DIR}
EXONERATE_OUTPUT_FILE="${OUT_DIR}exonerate_ferret_ncbi_pep_${SLURM_ARRAY_TASK_ID}.out"

exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}