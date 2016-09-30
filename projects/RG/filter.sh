#!/usr/bin/env bash

#for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123; do mkdir -p /home/genomerussia/main/fastq/filtered/GR00/${SAMPLE}; mkdir -p /home/genomerussia/main/analysis/stat/filtering/GR00/; mkdir -p /home/genomerussia/main/analysis/stat/filtering/GR00/${SAMPLE} ; FILES=($(ls /home/genomerussia/main/fastq/unpacked/GR00/${SAMPLE}/* | sed 's/.gz//')); echo ${FILES[0]}; echo ${FILES[1]}; /home/genomerussia/tools/Facut/bin/filter_by_mean_quality -t 20 -q phred33 -n illumina  -f ${FILES[0]} -r ${FILES[1]} -p /home/genomerussia/main/fastq/filtered/GR00/${SAMPLE}/${SAMPLE} > /home/genomerussia/main/analysis/stat/filtering/GR00/${SAMPLE}/${SAMPLE}.filtering.stat; done

for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123;
    do
    mkdir -p /home/genomerussia/main/fastq/filtered/GR00/${SAMPLE};
    mkdir -p /home/genomerussia/main/analysis/stat/filtering/GR00/;
    mkdir -p /home/genomerussia/main/analysis/stat/filtering/GR00/${SAMPLE};
    FILES=($(ls /home/genomerussia/main/fastq/unpacked/GR00/${SAMPLE}/* | sed 's/.gz//'));
    echo ${FILES[0]};
    echo ${FILES[1]};
    /home/genomerussia/tools/Facut/bin/filter_by_mean_quality -t 20 -q phred33 -n illumina  -f ${FILES[0]} -r ${FILES[1]} -p /home/genomerussia/main/fastq/filtered/GR00/${SAMPLE}/${SAMPLE} > /home/genomerussia/main/analysis/stat/filtering/GR00/${SAMPLE}/${SAMPLE}.filtering.stat;
    done
