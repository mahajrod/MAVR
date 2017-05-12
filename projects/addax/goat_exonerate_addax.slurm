#!/usr/bin/env bash
#SBATCH --array=1-640%64
#SBATCH --time=1000:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=8096       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=goat_addax
#SBATCH --error=/data4/genart/barefoot/addax/addax_annotation/err/goat_addax.%A_%a.err
#SBATCH --output=/data4/genart/barefoot/addax/addax_annotation/log/goat_addax.%A_%a.out

SPLITED_PEP_DIR="/data4/genart/barefoot/addax/addax_annotation/splited/"

EXONERATE_OUTPUT_DIR="/data4/genart/barefoot/addax/addax_annotation/output/"

GENOME="/data4/genart/barefoot/addax/addax_annotation/genome/K63g15d5e5R.scafSeq.gapclosed3.fas"

PEP_FILE="${SPLITED_PEP_DIR}/goat_pep_${SLURM_ARRAY_TASK_ID}.fasta"
EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/goat_exonerate_${SLURM_ARRAY_TASK_ID}.out"

exonerate --model protein2genome --showalignment --showquerygff --showtargetgff \
          -n 5 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}
