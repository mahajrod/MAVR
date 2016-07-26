#!/usr/bin/env bash
#SBATCH --array=117,166,175,220,304,313,358,389,436,469,488,496,535,572,575,601,652,67,682,692,722,738,766,76,773,843,860,866,91,95
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=12096       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Cats_ferret
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/cat_comparison/logs/exonerate_ferret_felcat.ensembl.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/cat_comparison/logs/exonerate_ferret_felcat.ensembl.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

REFERENCE_SPECIES="ferret"
SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/ferret/pep/ensembl/splited_pep/"
PEP_FILE="${SPLITED_PEP_DIR}mustela_putoris_furo_ensembl_${SLURM_ARRAY_TASK_ID}.fasta"

CAT_GENOMES_DIR="/work/pavlov/okochenova/genomes/cat_project/"
WORK_DIR="/work/pavlov/okochenova/genomes_workdir/cat_comparison/"

for SPECIES in FelCat_8.0; # AciJub_1.0 CroCro_1.0 LynPar_1.0 PanPar_1.0;
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






