#!/usr/bin/env bash

cd ~/data/genomes/saccharomyces_cerevisiae/LAN210/LAN210_v0.10m

samtools view -q 20 -b  1279_SMY732-4-7.bowtie2.sorted.rmdup.read_groups.bam > 1279_SMY732-4-7.bowtie2.sorted.rmdup.read_groups.q20.bam
samtools index 1279_SMY732-4-7.bowtie2.sorted.rmdup.read_groups.q20.bam

samtools faidx LAN210_v0.10m.fasta

java -jar ~/soft/PiCARD/picard.jar CreateSequenceDictionary R= LAN210_v0.10m.fasta O= LAN210_v0.10m.dict

cd ~/workdir/yeast/aksenova/illumina/alignment/1279_SMY732-4-7

~/workdir/yeast/aksenova/illumina/alignment/1279_SMY732-4-7$ java -jar ~/soft/PiCARD/picard.jar AddOrReplaceReadGroups I= 1279_SMY732-4-7.bowtie2.sorted.rmdup.bam O= 1279_SMY732-4-7.bowtie2.sorted.rmdup.read_groups.bam RGID=1 RGLB=1279_SMY732-4-7 RGPL=illumina RGPU=1 RGSM=1279_SMY732-4-7

samtools index 1279_SMY732-4-7.bowtie2.sorted.rmdup.read_groups.bam

cd ~/workdir/yeast/aksenova/illumina/SNPcall/1279_SMY732-4-7

java -jar ~/soft/GATK/GenomeAnalysisTK.jar -R ~/data/genomes/saccharomyces_cerevisiae/LAN210/LAN210_v0.10m/LAN210_v0.10m.fasta -T HaplotypeCaller  -I ../../alignment/1279_SMY732-4-7/1279_SMY732-4-7.bowtie2.sorted.rmdup.read_groups.bam -o 1279_SMY732-4-7.bowtie2.sorted.rmdup.raw.vcf -stand_emit_conf 20 -stand_call_conf 40 -nct 10 --sample_ploidy 1 -fixMisencodedQuals



cd ~/workdir/yeast/aksenova/torrent/alignment

~/soft/bwa-0.7.12/bwa mem -t 10 ~/data/genomes/saccharomyces_cerevisiae/LAN210/LAN210_v0.10m/bwa_index/LAN210_v0.10m ../fastq/filtered/01-2016_06_08_Aksenova.trimmed.quality_filtered.se.fq > 01-2016_06_08_Aksenova.trimmed.quality_filtered.se.bam &
~/soft/bwa-0.7.12/bwa mem -t 10 ~/data/genomes/saccharomyces_cerevisiae/LAN210/LAN210_v0.10m/bwa_index/LAN210_v0.10m ../fastq/filtered/02-2016_06_08_Aksenova.trimmed.quality_filtered.se.fq > 02-2016_06_08_Aksenova.trimmed.quality_filtered.se.bam &

tmap mapall -n 10 -f ~/data/genomes/saccharomyces_cerevisiae/LAN210/LAN210_v0.10m/LAN210_v0.10m.fasta -i fq -r ../fastq/filtered/01-2016_06_08_Aksenova.trimmed.quality_filtered.se.fq  -s 01-2016_06_08_Aksenova.trimmed.quality_filtered.se.tmap.bam -o 2 -v stage1 map1 map2 map3

tmap mapall -n 10 -f ~/data/genomes/saccharomyces_cerevisiae/LAN210/LAN210_v0.10m/LAN210_v0.10m.fasta -i fq -r ../fastq/filtered/02-2016_06_08_Aksenova.trimmed.quality_filtered.se.fq  -s 02-2016_06_08_Aksenova.trimmed.quality_filtered.se.tmap.bam -o 2 -v stage1 map1 map2 map3

for BAM in *.bam; do samtools sort -o ${BAM%%.*}.sorted.bam ${BAM}; samtools rmdup ${BAM%%.*}.sorted.bam -o ${BAM%%.*}.sorted.rmdup.bam;  done

for BAM in *.bam; do GR=${BAM%%.*}; PR=`basename ${BAM} .bam`; samtools sort ${BAM} | samtools rmdup - ${PR}.sorted.rmdup.bam; java -jar ~/soft/PiCARD/picard.jar AddOrReplaceReadGroups I= ${PR}.sorted.rmdup.bam O= ${PR}.sorted.rmdup.read_groups.bam RGID=1 RGLB=${GR} RGPL=iontorrent RGPU=1 RGSM=${GR}; samtools index ${PR}.sorted.rmdup.read_groups.bam; done
for BAM in *.read_groups.bam; do PR=`basename ${BAM} .read_groups.bam`; java -jar ~/soft/GATK/GenomeAnalysisTK.jar -R ~/data/genomes/saccharomyces_cerevisiae/LAN210/LAN210_v0.10m/LAN210_v0.10m.fasta -T HaplotypeCaller  -I ${BAM} -o ../SNPcall/${PR}.raw.vcf -stand_emit_conf 20 -stand_call_conf 40 -nct 10 --sample_ploidy 1; done
