#!/usr/bin/env bash
#SBATCH --array=1-920%920
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Dval_python
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/Dval/logs/exonerate_python.ncbi.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/Dval/logs/exonerate_python.ncbi.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/python/pep/python_ncbi/"
PEP_FILE="${SPLITED_PEP_DIR}python_ncbi_${SLURM_ARRAY_TASK_ID}.fasta"

GENOME="/work/pavlov/okochenova/genomes/Darevskia/D_val_v1.fa"

EXONERATE_OUTPUT_DIR="/work/pavlov/okochenova/genomes_workdir/Dval/exonerate_python_pep/"
mkdir -p ${EXONERATE_OUTPUT_DIR}
EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/exonerate_python_pep_${SLURM_ARRAY_TASK_ID}.out"

exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}

