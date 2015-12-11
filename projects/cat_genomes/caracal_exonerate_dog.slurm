#!/usr/bin/env bash
#SBATCH --array=1-992%350
#SBATCH --time=100:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=8192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=cara_dog
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/caracal/logs/exonerate_dog.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/caracal/logs/exonerate_dog.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/dog/pep/splited_pep/"

PEP_FILE="${SPLITED_PEP_DIR}canis_familiaris.sorted.description.uniq.no_selenocystein_${SLURM_ARRAY_TASK_ID}.fasta"
CARACAL_GENOME="/work/pavlov/okochenova/genomes/caracal/final.assembly.fasta"

mkdir -p /work/pavlov/okochenova/genomes_workdir/caracal/exonerate_dog_pep/
EXONERATE_OUTPUT_FILE="/work/pavlov/okochenova/genomes_workdir/caracal/exonerate_dog_pep/exonerate_dog_pep_${SLURM_ARRAY_TASK_ID}.out"

exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${CARACAL_GENOME} > ${EXONERATE_OUTPUT_FILE}