#!/usr/bin/env bash
#SBATCH --array=1-941%941
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Drad_gekko
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/Drad/logs/exonerate_gekko.ncbi.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/Drad/logs/exonerate_gekko.ncbi.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/gekko/pep/gekko_ncbi/"
PEP_FILE="${SPLITED_PEP_DIR}geckko_ncbi_${SLURM_ARRAY_TASK_ID}.fasta"

GENOME="/work/pavlov/okochenova/genomes/Darevskia/D_rad_v1.fa"

EXONERATE_OUTPUT_DIR="/work/pavlov/okochenova/genomes_workdir/Drad/exonerate_gekko_pep/"
mkdir -p ${EXONERATE_OUTPUT_DIR}
EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/exonerate_gekko_pep_${SLURM_ARRAY_TASK_ID}.out"

exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}

