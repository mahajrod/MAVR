#!/usr/bin/env bash

~/soft/MAVR/scripts/filter/filtering_mirna_pipeline.py -d fastq -o ./ -q 15 -f ~/soft/FaCut/bin/ -a ~/data/service_seq/trueseq_adapters_with_rev_com.fasta

for SAMPLE in `ls fastq`; do ~/soft/MAVR/scripts/filter/filtering_mirna_pipeline.py -d fastq -o ./ -q 15 -f ~/soft/FaCut/bin/ -a ~/data/service_seq/trueseq_adapters_with_rev_com.fasta -e 30 -s ${SAMPLE} & done

~/soft/MAVR/scripts/transcriptome/align_miRNA_by_star.py  -d filtered/final/ -o star_miRNA_hg_2mismatches_alignment/ -g reference/human_genome/ -t 50   -r ~/soft/STAR/bin/Linux_x86_64/ -z 0.1










for SAMPLE in `ls filtered/final/`; do mkdir -p alignment/mirBase.hsa/${SAMPLE}; bowtie2 -p 32 -x /home/skliver//workdir/preclampsia/reference/mirBase/mature.hsa.dna  --very-sensitive-local -U filtered/final/${SAMPLE}/${SAMPLE}.se.fq | samtools view -b  -@ 10 > alignment/mirBase.hsa/${SAMPLE}/${SAMPLE}.bam; done


for SAMPLE in `ls filtered/final/`; do mkdir -p alignment/mirBase.hsa/${SAMPLE}; ~/soft//bwa-0.7.12/bwa aln -o 0 -t 32 reference/mirBase/mature.hsa.dna.fa  filtered/final/${SAMPLE}/${SAMPLE}.se.fq | samtools view -b  -@ 10 > alignment/mirBase.hsa/${SAMPLE}/${SAMPLE}.bwa.bam 2>alignment/mirBase.hsa/${SAMPLE}/${SAMPLE}.alignment.bwa.log; done


