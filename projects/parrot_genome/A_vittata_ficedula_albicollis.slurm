#!/usr/bin/env bash
#SBATCH --array=1-942%942
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Avit_ficedula_albicollis
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/Avit/logs/exonerate_ficedula_albicollis.ensembl.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/Avit/logs/exonerate_ficedula_albicollis.ensembl.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/parrot/Ficedula_albicollis/splited/"
PEP_FILE="${SPLITED_PEP_DIR}Ficedula_albicollis_${SLURM_ARRAY_TASK_ID}.fasta"

GENOME="/work/pavlov/okochenova/genomes/parrot/Amazona_vittata/L_RNA_scaffolder.fasta"

EXONERATE_OUTPUT_DIR="/work/pavlov/okochenova/genomes_workdir/Avit/exonerate_ficedula_albicollisn_pep/"
mkdir -p ${EXONERATE_OUTPUT_DIR}
EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/exonerate_ficedula_albicollis_pep_${SLURM_ARRAY_TASK_ID}.out"

exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}
