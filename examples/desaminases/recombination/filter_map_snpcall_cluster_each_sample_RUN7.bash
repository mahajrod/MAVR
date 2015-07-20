#!/usr/bin/env bash
#SBATCH --array=1-38%5
#SBATCH -n 10
#SBATCH --time=100:00:00          # Run time in hh:mm:ss
#SBATCH --mem-per-cpu=8096       # Minimum memory required per CPU (in megabytes)
#SBATCH --job-name=filtering_ns
#SBATCH --error=/work/pavlov/okochenova/job_reports/RUN7/filtering_ns/job.%A_%a.err
#SBATCH --output=/work/pavlov/okochenova/job_reports/RUN7/filtering_ns/job.%A_%a.out 

module load compiler/gcc/4.8 fastx_toolkit/0.0.14 samtools/0.1 bowtie/2.2 python/2.7 cutadapt/1.4

source /work/pavlov/okochenova/profile
SAMPLES=(N085-LAN210-Can-PmCDA1-NA-RUN7-D1 N086-LAN210-Can-PmCDA1-NA-RUN7-D1 N087-LAN210-Can-PmCDA1-NA-RUN7-D1 N088-LAN210-Can-PmCDA1-NA-RUN7-D1 N089-LAN210-Can-PmCDA1-NA-RUN7-D1 N090-LAN210-Can-PmCDA1-NA-RUN7-D1 N091-LAN210-Can-PmCDA1-NA-RUN7-D3 N092-LAN210-Can-PmCDA1-NA-RUN7-D6 N093-LAN210-Can-PmCDA1-NA-RUN7-D6 N094-LAN210-Can-PmCDA1-NA-RUN7-D6 N095-LAN210-Can-PmCDA1-NA-RUN7-D6 N096-LAN210-Can-PmCDA1-NA-RUN7-D6 N097-LAN210-Can-AID-NA-RUN7-D1 N098-LAN210-Can-AID-NA-RUN7-D1 N100-LAN210-Can-AID-NA-RUN7-D3 N101-LAN210-Can-AID-NA-RUN7-D3 N102-LAN210-Can-AID-NA-RUN7-D3 N103-LAN210-Can-AID-NA-RUN7-D6 N104-LAN210-Can-AID-NA-RUN7-D6 N105-LAN210-Can-AID-NA-RUN7-D6 N106-LAN210-Can-AID-NA-RUN7-D6 N107-LAN210-Can-A1-NA-RUN7-D1 N108-LAN210-Can-A1-NA-RUN7-D1 N109-LAN210-Can-A1-NA-RUN7-D1 N110-LAN210-Can-A1-NA-RUN7-D1 N111-LAN210-Can-A1-NA-RUN7-D1 N112-LAN210-Can-A1-NA-RUN7-D1 N113-LAN210-Can-A1-NA-RUN7-D3 N114-LAN210-Can-A1-NA-RUN7-D3 N115-LAN210-Can-A1-NA-RUN7-D3 N116-LAN210-Can-A1-NA-RUN7-D3 N117-LAN210-Can-A1-NA-RUN7-D3 N118-LAN210-Can-A1-NA-RUN7-D3 N120-LAN210-Can-A1-NA-RUN7-D6 N121-LAN210-Can-A1-NA-RUN7-D6 N122-LAN210-Can-A1-NA-RUN7-D6 N123-LAN210-Can-A1-NA-RUN7-D6 N124-LAN210-Can-A1-NA-RUN7-D6)

cd /work/pavlov/okochenova/fastq/RUN7/
PICARD_DIR="/work/pavlov/okochenova/soft/picard-tools-1.115/picard-tools-1.115"
let "SAMPLE_INDEX=${SLURM_ARRAY_TASK_ID}-1"
CURRENT_SAMPLE=${SAMPLES[${SAMPLE_INDEX}]} 
cd ${CURRENT_SAMPLE}
echo ${CURRENT_SAMPLE}

mkdir -p filtered_ns
cd filtered_ns
/work/pavlov/okochenova/soft/MAVR/scripts/filter/remove_terminal_ns.py -l ../${CURRENT_SAMPLE}_1.fastq -r ../${CURRENT_SAMPLE}_2.fastq -o ./ -m 30

fastqc --nogroup -t 4 *.f*q

mkdir -p filtered

/work/pavlov/okochenova/soft/rm_reads --adapters /work/pavlov/okochenova/kmers/service/illumina_testseq_trueseq_with_rc_23_mer.kmer -1 ${CURRENT_SAMPLE}_1.filtered_1.fastq -2 ${CURRENT_SAMPLE}_2.filtered_2.fastq -o ./filtered

cd filtered

fastqc --nogroup -t 4 *.f*q

mkdir -p trimmed

trim_galore --paired --length 30 --phred33 -q 20 -o ./trimmed --dont_gzip -a CTGTCTCTTATACACATCT ${CURRENT_SAMPLE}_1.ok.fastq ${CURRENT_SAMPLE}_2.ok.fastq
trim_galore --length 30 --phred33 -q 20 -o ./trimmed --dont_gzip -a CTGTCTCTTATACACATCT ${CURRENT_SAMPLE}_1.se.fastq
trim_galore --length 30 --phred33 -q 20 -o ./trimmed --dont_gzip -a CTGTCTCTTATACACATCT ${CURRENT_SAMPLE}_2.se.fastq

cd trimmed

fastqc -t 4 --nogroup *.f*q

/work/pavlov/okochenova/soft/MAVR/scripts/alignment/map_reads.py -i /work/pavlov/okochenova/reference/LAN210/LAN210_v0.10m/LAN210_v0.10m -a bowtie2 -t 10 -b /work/pavlov/okochenova/reference/LAN210/LAN210_v0.10m/chromosomes.bed -r ${CURRENT_SAMPLE}_1.ok_val_1.fq -l ${CURRENT_SAMPLE}_2.ok_val_2.fq -u ${CURRENT_SAMPLE}_1.se_trimmed.fq,${CURRENT_SAMPLE}_2.se_trimmed.fq -q phred33 -n  -d ${PICARD_DIR} -p ${CURRENT_SAMPLE}

ALIGNMENT_FILE=${CURRENT_SAMPLE}_final_with_groups.bam
/work/pavlov/okochenova/soft/MAVR/scripts/snpcall/snpcall.py -b ${ALIGNMENT_FILE} -p ${CURRENT_SAMPLE} -r /work/pavlov/okochenova/reference/LAN210/LAN210_v0.10m/LAN210_v0.10m.fasta -t 10 -q 40 -u 100 -g /work/pavlov/okochenova/soft/GATK/

/work/pavlov/okochenova/soft/MACE/scripts/clustering_pipeline.py -r /work/pavlov/okochenova/reference/LAN210/LAN210_v0.10m/LAN210_v0.10m.fasta -a /work/pavlov/okochenova/reference/LAN210/LAN210_v0.10m/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3  -m /work/pavlov/okochenova/reference/LAN210/LAN210_v0.10m/LAN210_v0.10m_masked_all_not_in_good_genes.gff -s ${CURRENT_SAMPLE}_GATK_best_SNP
