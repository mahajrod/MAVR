#!/usr/bin/env bash

~/soft/MAVR/scripts/filter/filtering_mirna_pipeline.py -d fastq -o ./ -q 15 -f ~/soft/FaCut/bin/ -a ~/data/service_seq/trueseq_adapters_with_rev_com.fasta

for SAMPLE in `ls fastq`; do ~/soft/MAVR/scripts/filter/filtering_mirna_pipeline.py -d fastq -o ./ -q 15 -f ~/soft/FaCut/bin/ -a ~/data/service_seq/trueseq_adapters_with_rev_com.fasta -e 30 -s ${SAMPLE} & done





~/soft/MAVR/scripts/transcriptome/align_miRNA_by_star.py  -d filtered/final/ -o star_miRNA_hg_2mismatches_alignment/ -g reference/human_genome/ -t 50   -r ~/soft/STAR/bin/Linux_x86_64/ -z 0.1




#skliver@supermicro:~/workdir/preclampsia/star_miRNA_hg_2mismatches_alignment$
for SP in `ls`; do ~/soft/MAVR/scripts/transcriptome/miRNA/count_miRNA_by_subread.py -a ${SP}/Aligned.sortedByCoord.out.bam -e ../reference/human_genome/GCF_000001405.26_GRCh38_genomic.renamed.ncrna.with_exons.only_miRNA.gff -t 1 --feature_id_attribute product -s ${SP} -o ${SP}/${SP} & done

~/soft/MAVR/scripts/transcriptome/differential_expression/combine_count_file.py -f 100A/100A.all_adjusted_reads.count,113A/113A.all_adjusted_reads.count,115A/115A.all_adjusted_reads.count,115p/115p.all_adjusted_reads.count,117A/117A.all_adjusted_reads.count,124p/124p.all_adjusted_reads.count,136p/136p.all_adjusted_reads.count,139p/139p.all_adjusted_reads.count,147p/147p.all_adjusted_reads.count,163p/163p.all_adjusted_reads.count,169p/169p.all_adjusted_reads.count,185p/185p.all_adjusted_reads.count,187p/187p.all_adjusted_reads.count,188p/188p.all_adjusted_reads.count,20p/20p.all_adjusted_reads.count,33A/33A.all_adjusted_reads.count,73A/73A.all_adjusted_reads.count,74p/74p.all_adjusted_reads.count,78A/78A.all_adjusted_reads.count,84A/84A.all_adjusted_reads.count,86A/86A.all_adjusted_reads.count,88A/88A.all_adjusted_reads.count,92A/92A.all_adjusted_reads.count,97A/97A.all_adjusted_reads.count -s 100A,113A,115A,115p,117A,124p,136p,139p,147p,163p,169p,185p,187p,188p,20p,33A,73A,74p,78A,84A,86A,88A,92A,97A -o all.samples.counts

for SP in `ls`; do ~/soft/MAVR/scripts/transcriptome/miRNA/count_miRNA_by_subread.py -a ${SP}/Aligned.sortedByCoord.out.bam -e ../reference/human_genome/GCF_000001405.26_GRCh38_genomic.renamed.ncrna.with_exons.only_miRNA.gff -t 1 -m 0.95 --feature_id_attribute product -s ${SP} -o ${SP}/${SP}.overlap_95 & done

~/soft/MAVR/scripts/transcriptome/differential_expression/combine_count_file.py -f 100A/100A.overlap_95.all_adjusted_reads.count,113A/113A.overlap_95.all_adjusted_reads.count,115A/115A.overlap_95.all_adjusted_reads.count,115p/115p.overlap_95.all_adjusted_reads.count,117A/117A.overlap_95.all_adjusted_reads.count,124p/124p.overlap_95.all_adjusted_reads.count,136p/136p.overlap_95.all_adjusted_reads.count,139p/139p.overlap_95.all_adjusted_reads.count,147p/147p.overlap_95.all_adjusted_reads.count,163p/163p.overlap_95.all_adjusted_reads.count,169p/169p.overlap_95.all_adjusted_reads.count,185p/185p.overlap_95.all_adjusted_reads.count,187p/187p.overlap_95.all_adjusted_reads.count,188p/188p.overlap_95.all_adjusted_reads.count,20p/20p.overlap_95.all_adjusted_reads.count,33A/33A.overlap_95.all_adjusted_reads.count,73A/73A.overlap_95.all_adjusted_reads.count,74p/74p.overlap_95.all_adjusted_reads.count,78A/78A.overlap_95.all_adjusted_reads.count,84A/84A.overlap_95.all_adjusted_reads.count,86A/86A.overlap_95.all_adjusted_reads.count,88A/88A.overlap_95.all_adjusted_reads.count,92A/92A.overlap_95.all_adjusted_reads.count,97A/97A.overlap_95.all_adjusted_reads.count -s 100A,113A,115A,115p,117A,124p,136p,139p,147p,163p,169p,185p,187p,188p,20p,33A,73A,74p,78A,84A,86A,88A,92A,97A -o all.samples.overlap_95.counts


for SP in `ls`; do ~/Genetics/MAVR/scripts/transcriptome/miRNA/count_miRNA_by_subread.py -a ${SP}/Aligned.sortedByCoord.out.bam -e /home/mahajrod/Genetics/ksi/preclampsia/GCF_000001405.26_GRCh38_genomic.renamed.ncrna.with_exons.only_miRNA_added_mirBase.gff -t 3 -m 0.90 --feature_id_attribute product -s ${SP} -o ${SP}/${SP}.overlap_90 ; done


~/Genetics/MAVR/scripts/transcriptome/differential_expression/combine_count_file.py -f 100A/100A.overlap_95.all_adjusted_reads.count,113A/113A.overlap_95.all_adjusted_reads.count,115A/115A.overlap_95.all_adjusted_reads.count,115p/115p.overlap_95.all_adjusted_reads.count,117A/117A.overlap_95.all_adjusted_reads.count,124p/124p.overlap_95.all_adjusted_reads.count,136p/136p.overlap_95.all_adjusted_reads.count,139p/139p.overlap_95.all_adjusted_reads.count,147p/147p.overlap_95.all_adjusted_reads.count,163p/163p.overlap_95.all_adjusted_reads.count,169p/169p.overlap_95.all_adjusted_reads.count,185p/185p.overlap_95.all_adjusted_reads.count,187p/187p.overlap_95.all_adjusted_reads.count,188p/188p.overlap_95.all_adjusted_reads.count,20p/20p.overlap_95.all_adjusted_reads.count,33A/33A.overlap_95.all_adjusted_reads.count,73A/73A.overlap_95.all_adjusted_reads.count,74p/74p.overlap_95.all_adjusted_reads.count,78A/78A.overlap_95.all_adjusted_reads.count,84A/84A.overlap_95.all_adjusted_reads.count,86A/86A.overlap_95.all_adjusted_reads.count,88A/88A.overlap_95.all_adjusted_reads.count,92A/92A.overlap_95.all_adjusted_reads.count,97A/97A.overlap_95.all_adjusted_reads.count -s 100A,113A,115A,115p,117A,124p,136p,139p,147p,163p,169p,185p,187p,188p,20p,33A,73A,74p,78A,84A,86A,88A,92A,97A -o all.samples.overlap_95.counts

#using miRbase annotations
for SP in `ls`; do ~/Genetics/MAVR/scripts/transcriptome/miRNA/count_miRNA_by_subread.py -a ${SP}/Aligned.sortedByCoord.out.bam -e /home/mahajrod/Databases/Mirbase/hsa.ncbi.names.gff3 -t 3 -m 0.90 --feature_type miRNA --feature_id_attribute Name -s ${SP} -o ${SP}/${SP}.mirbase.overlap_90 ; done

~/Genetics/MAVR/scripts/transcriptome/differential_expression/combine_count_file.py -f 100A/100A.mirbase.overlap_90.all_adjusted_reads.count,113A/113A.mirbase.overlap_90.all_adjusted_reads.count,115A/115A.mirbase.overlap_90.all_adjusted_reads.count,115p/115p.mirbase.overlap_90.all_adjusted_reads.count,117A/117A.mirbase.overlap_90.all_adjusted_reads.count,124p/124p.mirbase.overlap_90.all_adjusted_reads.count,136p/136p.mirbase.overlap_90.all_adjusted_reads.count,139p/139p.mirbase.overlap_90.all_adjusted_reads.count,147p/147p.mirbase.overlap_90.all_adjusted_reads.count,163p/163p.mirbase.overlap_90.all_adjusted_reads.count,169p/169p.mirbase.overlap_90.all_adjusted_reads.count,185p/185p.mirbase.overlap_90.all_adjusted_reads.count,187p/187p.mirbase.overlap_90.all_adjusted_reads.count,188p/188p.mirbase.overlap_90.all_adjusted_reads.count,20p/20p.mirbase.overlap_90.all_adjusted_reads.count,33A/33A.mirbase.overlap_90.all_adjusted_reads.count,73A/73A.mirbase.overlap_90.all_adjusted_reads.count,74p/74p.mirbase.overlap_90.all_adjusted_reads.count,78A/78A.mirbase.overlap_90.all_adjusted_reads.count,84A/84A.mirbase.overlap_90.all_adjusted_reads.count,86A/86A.mirbase.overlap_90.all_adjusted_reads.count,88A/88A.mirbase.overlap_90.all_adjusted_reads.count,92A/92A.mirbase.overlap_90.all_adjusted_reads.count,97A/97A.mirbase.overlap_90.all_adjusted_reads.count -s 100A,113A,115A,115p,117A,124p,136p,139p,147p,163p,169p,185p,187p,188p,20p,33A,73A,74p,78A,84A,86A,88A,92A,97A -o all.samples.mirbase.overlap_90.counts



for SAMPLE in `ls filtered/final/`; do mkdir -p alignment/mirBase.hsa/${SAMPLE}; bowtie2 -p 32 -x /home/skliver//workdir/preclampsia/reference/mirBase/mature.hsa.dna  --very-senonly_miRNAsitive-local -U filtered/final/${SAMPLE}/${SAMPLE}.se.fq | samtools view -b  -@ 10 > alignment/mirBase.hsa/${SAMPLE}/${SAMPLE}.bam; done


for SAMPLE in `ls filtered/final/`; do mkdir -p alignment/mirBase.hsa/${SAMPLE}; ~/soft//bwa-0.7.12/bwa aln -o 0 -t 32 reference/mirBase/mature.hsa.dna.fa  filtered/final/${SAMPLE}/${SAMPLE}.se.fq | samtools view -b  -@ 10 > alignment/mirBase.hsa/${SAMPLE}/${SAMPLE}.bwa.bam 2>alignment/mirBase.hsa/${SAMPLE}/${SAMPLE}.alignment.bwa.log; done


