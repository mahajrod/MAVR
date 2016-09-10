#!/usr/bin/env bash

filtration


cd /mnt/guatemala/skliver/Amazona_vittata/transcriptome/filtered_trimmomatic



#third variant
cd /mnt/guatemala/skliver/Amazona_vittata/transcriptome/filtered_trimmomatic_only
for SAMPLE in MAD-64-05 MAD-64-06 MAD-64-07 PLI-59-T0853; do FILE=`ls /mnt/guatemala/skliver/Amazona_vittata/transcriptome/se/*${SAMPLE}*`; java -jar /storage1/home//skliver/Soft/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 10 ${FILE} ${SAMPLE}.trimmomatic.slw.10.15.ml.50.fq SLIDINGWINDOW:10:15 MINLEN:50 > ${SAMPLE}.trimmomatic.slw.10.15.ml.50.log 2>&1; done

# create index
STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ../L_RNA_scaffolder.fasta  --runThreadN 25

#alignment_1_pass
cd /mnt/guatemala/skliver/Amazona_vittata/transcriptome/alignment_trimmomatic_only/

for SAMPLE in MAD-64-05 MAD-64-06 MAD-64-07 PLI-59-T0853; do mkdir -p ${SAMPLE}; cd ${SAMPLE}; echo "Handling ${SAMPLE}"; ~/Soft/STAR/bin/Linux_x86_64/STAR --runThreadN 20 --genomeDir /mnt/guatemala/skliver/Amazona_vittata/assembly/star_only_genome_idx/ --readFilesIn /mnt/guatemala/skliver/Amazona_vittata/transcriptome/filtered_trimmomatic_only/${SAMPLE}.trimmomatic.slw.10.15.ml.50.fq --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic;  cd ../; done

#alignmenet_2_pass all samples simalteniously

 ~/Soft/STAR/bin/Linux_x86_64/STAR --runThreadN 20 --genomeDir /mnt/guatemala/skliver/Amazona_vittata/assembly/star_only_genome_idx/ --readFilesIn /mnt/guatemala/skliver/Amazona_vittata/transcriptome/filtered_trimmomatic_only/MAD-64-06.trimmomatic.slw.10.15.ml.50.fq,/mnt/guatemala/skliver/Amazona_vittata/transcriptome/filtered_trimmomatic_only/MAD-64-07.trimmomatic.slw.10.15.ml.50.fq,/mnt/guatemala/skliver/Amazona_vittata/transcriptome/filtered_trimmomatic_only/MAD-64-05.trimmomatic.slw.10.15.ml.50.fq,/mnt/guatemala/skliver/Amazona_vittata/transcriptome/filtered_trimmomatic_only/PLI-59-T0853.trimmomatic.slw.10.15.ml.50.fq --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --twopassMode Basic --sjdbFileChrStartEnd /mnt/guatemala/skliver/Amazona_vittata/transcriptome/alignment_trimmomatic_only/MAD-64-05/_STARpass1/SJ.out.tab /mnt/guatemala/skliver/Amazona_vittata/transcriptome/alignment_trimmomatic_only/MAD-64-06/_STARpass1/SJ.out.tab /mnt/guatemala/skliver/Amazona_vittata/transcriptome/alignment_trimmomatic_only/MAD-64-07/_STARpass1/SJ.out.tab /mnt/guatemala/skliver/Amazona_vittata/transcriptome/alignment_trimmomatic_only/PLI-59-T0853/_STARpass1/SJ.out.tab
