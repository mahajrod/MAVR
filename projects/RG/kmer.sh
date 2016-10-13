#!/usr/bin/env bash

#for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123; do mkdir -p  /home/genomerussia/main/analysis/jf/GR00/${SAMPLE}; PYTHONPATH=${PYTHONPATH}:/home/genomerussia/tools/MAVR/ /home/genomerussia/tools/MAVR/scripts/kmer/draw_kmer_distribution_from_fastq.py -m 23 -t 60 -b -s 30G -e png -o /home/genomerussia/main/analysis/jf/GR00/${SAMPLE}/${SAMPLE} -w 3 -g 80 -i `ls -m /home/genomerussia/main/fastq/unpacked/GR00/${SAMPLE}/* | sed -r "s/, /,/g" | tr -d '\n'`;  done

for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123;
    do
    mkdir -p /home/genomerussia/main/analysis/jf/GR00/${SAMPLE};
    PYTHONPATH=${PYTHONPATH}:/home/genomerussia/tools/MAVR/ /home/genomerussia/tools/MAVR/scripts/kmer/draw_kmer_distribution_from_fastq.py -m 23 -t 60 -b -s 30G -e png -o /home/genomerussia/main/analysis/jf/GR00/${SAMPLE}/${SAMPLE} -w 3 -g 80 -i `ls -m /home/genomerussia/main/fastq/unpacked/GR00/${SAMPLE}/* | sed -r "s/, /,/g" | tr -d '\n'`;
    done

#for SAMPLE in GR0082 GR0081 GR0122 GR0121 GR0123; do mkdir -p /home/skliver/tmp/kmer/GR00/${SAMPLE}/; FILES=($(ls /home/skliver/tmp/GR00/${SAMPLE}/* | sed 's/.gz//')); echo ${FILES[0]}; echo ${FILES[1]}; echo ${FILES[2]}; echo ${FILES[3]}; /home/genomerussia/tools/MAVR/scripts/kmer/draw_kmer_distribution_from_fastq.py -i ${FILES[0]},${FILES[1]} -o /home/skliver/tmp/kmer/GR00/${SAMPLE}/`basename ${FILES[0]} .fq` -m 23 -s 20G -t 30 -b -w 3 -g 80 -e png; /home/genomerussia/tools/MAVR/scripts/kmer/draw_kmer_distribution_from_fastq.py -i ${FILES[2]},${FILES[3]} -o /home/skliver/tmp/kmer/GR00/${SAMPLE}/`basename ${FILES[2]} .fq` -m 23 -s 20G -t 30 -b -w 3 -g 80 -e png; done