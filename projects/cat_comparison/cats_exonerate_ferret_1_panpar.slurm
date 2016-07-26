#!/usr/bin/env bash
#SBATCH --array=117,175,692,722,76
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=18096       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Cats_ferret
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/cat_comparison/logs/exonerate_ferret.ensembl.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/cat_comparison/logs/exonerate_ferret.ensembl.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

REFERENCE_SPECIES="ferret"
SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/ferret/pep/ensembl/splited_pep/"
PEP_FILE="${SPLITED_PEP_DIR}mustela_putoris_furo_ensembl_${SLURM_ARRAY_TASK_ID}.fasta"

CAT_GENOMES_DIR="/work/pavlov/okochenova/genomes/cat_project/"
WORK_DIR="/work/pavlov/okochenova/genomes_workdir/cat_comparison/"

for SPECIES in PanPar_1.0; # PanLeo_1.0 PanTig_1.0 PanOnc_1.0 FelCat_8.0;
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






