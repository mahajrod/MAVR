#!/usr/bin/env bash
#SBATCH --array=1-9%9
#SBATCH -n 2
#SBATCH --time=100:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=4096       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=clustering_sets
#SBATCH --error=/work/pavlov/okochenova/job_reports/RUN7/clustering_sets.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/job_reports/RUN7/clustering_sets.%A_%a.out

module load compiler/gcc/4.8 python/2.7

source /work/pavlov/okochenova/profile

SCRIPT=/work/pavlov/okochenova/soft/MACE/scripts/clustering_pipeline.py
WORKDIR=/work/pavlov/okochenova/combined_sets/
VCF_DIR=${WORKDIR}raw/
REFERENCE=/work/pavlov/okochenova/reference/LAN210/LAN210_v0.10m/LAN210_v0.10m.fasta
REFERENCE_ANNOTATIONS=/work/pavlov/okochenova/reference/LAN210/LAN210_v0.10m/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3
REFERENCE_MASKING=/work/pavlov/okochenova/reference/LAN210/LAN210_v0.10m/LAN210_v0.10m_masked_all_not_in_good_genes.gff

SAMPLE_SETS=(A1_D1 A1_D3 A1_D6 AID_D1 AID_D3 AID_D6 PmCDA1_D1 PmCDA1_D3 PmCDA1_D6)

cd ${WORKDIR}

SAMPLE_SET_INDEX=
let "SAMPLE_SET_INDEX=${SLURM_ARRAY_TASK_ID}-1"
CURRENT_SAMPLE_SET=${SAMPLE_SETS[${SAMPLE_SET_INDEX}]}
cd ${CURRENT_SAMPLE_SET}
echo ${CURRENT_SAMPLE_SET}

${SCRIPT} -r ${REFERENCE} -a ${REFERENCE_ANNOTATIONS} -m ${REFERENCE_MASKING} -f ${VCF_DIR}${CURRENT_SAMPLE_SET}_raw.vcf -s ${CURRENT_SAMPLE_SET} -y ${CURRENT_SAMPLE_SET}
