#!/usr/bin/env bash
#SBATCH --array=1-967%967
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Avit_taeniopygia_guttata
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/Avit/logs/exonerate_taeniopygia_guttata.ensembl.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/Avit/logs/exonerate_taeniopygia_guttata.ensembl.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/parrot/Taeniopygia_guttata/splited/"
PEP_FILE="${SPLITED_PEP_DIR}Taeniopygia_guttata_${SLURM_ARRAY_TASK_ID}.fasta"

GENOME="/work/pavlov/okochenova/genomes/parrot/Amazona_vittata/L_RNA_scaffolder.fasta"

EXONERATE_OUTPUT_DIR="/work/pavlov/okochenova/genomes_workdir/Avit/exonerate_taeniopygia_guttata_pep/"
mkdir -p ${EXONERATE_OUTPUT_DIR}
EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/exonerate_taeniopygia_guttata_pep_${SLURM_ARRAY_TASK_ID}.out"

exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}
