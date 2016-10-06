#!/usr/bin/env bash

#for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123; do mkdir -p /home/genomerussia/main/analysis/stat/adapters//GR00/; mkdir -p /home/genomerussia/main/analysis/stat/adapters//GR00/${SAMPLE}; FILES=($(ls /home/genomerussia/main/fastq/unpacked/GR00/${SAMPLE}/* | sed 's/.gz//')); echo ${FILES[0]}; echo ${FILES[1]}; /home/genomerussia/tools/Cookiecutter/src/counter -1 ${FILES[0]} -2 ${FILES[1]} -o /home/genomerussia/main/analysis/stat/adapters/GR00/${SAMPLE}/ -f /home/genomerussia/tools/service_sequences/trueseq_adapters_with_rev_com_23_mer.kmer > /home/genomerussia/main/analysis/stat/adapters/GR00/${SAMPLE}${SAMPLE}.adapters.stat; done


for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123;
    do
    mkdir -p /home/genomerussia/main/analysis/stat/adapters//GR00/;
    mkdir -p /home/genomerussia/main/analysis/stat/adapters//GR00/${SAMPLE};
    FILES=($(ls /home/genomerussia/main/fastq/unpacked/GR00/${SAMPLE}/* | sed 's/.gz//'));
    echo ${FILES[0]};
    echo ${FILES[1]};
    /home/genomerussia/tools/Cookiecutter/src/counter -1 ${FILES[0]} -2 ${FILES[1]} -o /home/genomerussia/main/analysis/stat/adapters/GR00/${SAMPLE}/ -f /home/genomerussia/tools/service_sequences/trueseq_adapters_with_rev_com_23_mer.kmer > /home/genomerussia/main/analysis/stat/adapters/GR00/${SAMPLE}${SAMPLE}.adapters.stat;
    done


for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123; do mkdir -p /home/skliver/tmp/GR00/${SAMPLE}; FILES=($(ls /home/genomerussia/main/fastq/unpacked/GR00/${SAMPLE}/* | sed 's/.gz//')); echo ${FILES[0]}; echo ${FILES[1]}; /home/genomerussia/tools/Facut/bin/split_fastq_by_lane  -n illumina -f ${FILES[0]} -r ${FILES[1]} -p /home/skliver/tmp/GR00/${SAMPLE}/${SAMPLE} ; done
