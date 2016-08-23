#!/usr/bin/env bash
#SBATCH --array=1-940%940
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4096       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Cats_human
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/cat_comparison/logs/exonerate_human.ncbi.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/cat_comparison/logs/exonerate_human.ncbi.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

REFERENCE_SPECIES="human"
SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/human/ncbi/pep/splited/"
PEP_FILE="${SPLITED_PEP_DIR}homo_sapiens.ncbi_${SLURM_ARRAY_TASK_ID}.fasta"

CAT_GENOMES_DIR="/work/pavlov/okochenova/genomes/cat_project/"
WORK_DIR="/work/pavlov/okochenova/genomes_workdir/cat_comparison/"

for SPECIES in CarCar_1.0 PriViv_1.0 PriBen_1.0;
    do
    TARGET_GENOME="${CAT_GENOMES_DIR}${SPECIES}.fixed.fa"
    TARGET_GENOME_DIR="${WORK_DIR}${SPECIES}/"
    EXONERATE_OUTPUT_DIR="${TARGET_GENOME_DIR}${REFERENCE_SPECIES}/"
    EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/exonerate_${REFERENCE_SPECIES}_pep_${SLURM_ARRAY_TASK_ID}.out"
    mkdir -p ${TARGET_GENOME_DIR}
    mkdir -p ${EXONERATE_OUTPUT_DIR}
    exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${TARGET_GENOME} > ${EXONERATE_OUTPUT_FILE}
    done






