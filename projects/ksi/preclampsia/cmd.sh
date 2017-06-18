#!/usr/bin/env bash

for SAMPLE in `ls filtered/final/`; do mkdir -p alignment/mirBase.hsa/${SAMPLE}; bowtie2 -p 32 -x /home/skliver//workdir/preclampsia/reference/mirBase/mature.hsa.dna  --very-sensitive-local -U filtered/final/${SAMPLE}/${SAMPLE}.se.fq | samtools view -b  -@ 10 > alignment/mirBase.hsa/${SAMPLE}/${SAMPLE}.bam; done


for SAMPLE in `ls filtered/final/`; do mkdir -p alignment/mirBase.hsa/${SAMPLE}; ~/soft//bwa-0.7.12/bwa aln -o 0 -t 32 reference/mirBase/mature.hsa.dna.fa  filtered/final/${SAMPLE}/${SAMPLE}.se.fq | samtools view -b  -@ 10 > alignment/mirBase.hsa/${SAMPLE}/${SAMPLE}.bwa.bam 2>alignment/mirBase.hsa/${SAMPLE}/${SAMPLE}.alignment.bwa.log; done


