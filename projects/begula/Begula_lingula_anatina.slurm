#!/usr/bin/env bash
#SBATCH --array=1-948%948
#SBATCH --time=100:00:00         # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4192       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=Beg_lin
#SBATCH --error=/work/pavlov/okochenova/genomes_workdir/begula/logs/exonerate_lingula.ensembl.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/genomes_workdir/begula/logs/exonerate_lingula.ensembl.%A_%a.out

module load exonerate/2.2
source /work/pavlov/okochenova/profile

SPLITED_PEP_DIR="/work/pavlov/okochenova/genomes/lingula_anatina/splited/"
PEP_FILE="${SPLITED_PEP_DIR}lingula_anatina_${SLURM_ARRAY_TASK_ID}.fasta"

GENOME="/work/pavlov/okochenova/genomes/begula/begual.contigs.fasta"

EXONERATE_OUTPUT_DIR="/work/pavlov/okochenova/genomes_workdir/begula/exonerate_lingula_anatina_pep/"
mkdir -p ${EXONERATE_OUTPUT_DIR}
EXONERATE_OUTPUT_FILE="${EXONERATE_OUTPUT_DIR}/exonerate_lingula_anatina_pep_${SLURM_ARRAY_TASK_ID}.out"

exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME} > ${EXONERATE_OUTPUT_FILE}

GENOME_UN="/work/pavlov/okochenova/genomes/begula/begual.unassembled.fasta"

EXONERATE_UN_OUTPUT_DIR="/work/pavlov/okochenova/genomes_workdir/begula/exonerate_un_lingula_anatina_pep/"
mkdir -p ${EXONERATE_UN_OUTPUT_DIR}
EXONERATE_UN_OUTPUT_FILE="${EXONERATE_UN_OUTPUT_DIR}/exonerate_lingula_anatina_pep_${SLURM_ARRAY_TASK_ID}.out"

exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
           -n 10 -q ${PEP_FILE} -t ${GENOME_UN} > ${EXONERATE_UN_OUTPUT_FILE}