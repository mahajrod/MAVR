#!/usr/bin/env bash
#SBATCH --array=1-100
#SBATCH --mem-per-cpu=8096       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=c-mul
#SBATCH --error=/data4/genart/barefoot/caridina/caridina_annotation/err/C-mul.%A_%a.err
#SBATCH --output=/data4/genart/barefoot/caridina/caridina_annotation/log/C-mul.%A_%a.out

#/nfs/mfnstore-lin/export/mfn-genom-1/KLIVER/workdir/shrimp/CAP3/C-mul/100+/splited          ----->>>>>
SPLITED_PEP_DIR="/data4/genart/barefoot/caridina/caridina_annotation/C-mul/splited/"

EXONERATE_OUTPUT_DIR="/data4/genart/barefoot/caridina/caridina_annotation/C-mul/output/"

#/nfs/mfnstore-lin/export/mfn-genom-1/KLIVER/workdir/caridina_multidentata/genome/assemblies/PEP_scaffolder_one_link.fasta     ----->>>>>
GENOME="/data4/genart/barefoot/caridina/caridina_annotation/genome/PEP_scaffolder_one_link.fasta"

PEP_FILE="${SPLITED_PEP_DIR}/C-mul_${SLURM_ARRAY_TASK_ID}.fasta"
EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/C-mul_${SLURM_ARRAY_TASK_ID}.out"

exonerate --model protein2genome --showalignment --showquerygff --showtargetgff \
          -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}
