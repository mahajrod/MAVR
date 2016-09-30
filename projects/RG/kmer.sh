#!/usr/bin/env bash

#for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123; do mkdir -p  /home/genomerussia/main/analysis/jf/GR00/${SAMPLE}; PYTHONPATH=${PYTHONPATH}:/home/genomerussia/tools/MAVR/ /home/genomerussia/tools/MAVR/scripts/kmer/draw_kmer_distribution_from_fastq.py -m 23 -t 60 -b -s 30G -e png -o /home/genomerussia/main/analysis/jf/GR00/${SAMPLE}/${SAMPLE} -w 3 -g 80 -i `ls -m /home/genomerussia/main/fastq/unpacked/GR00/${SAMPLE}/* | sed -r "s/, /,/g" | tr -d '\n'`;  done

for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123;
    do
    mkdir -p /home/genomerussia/main/analysis/jf/GR00/${SAMPLE};
    PYTHONPATH=${PYTHONPATH}:/home/genomerussia/tools/MAVR/ /home/genomerussia/tools/MAVR/scripts/kmer/draw_kmer_distribution_from_fastq.py -m 23 -t 60 -b -s 30G -e png -o /home/genomerussia/main/analysis/jf/GR00/${SAMPLE}/${SAMPLE} -w 3 -g 80 -i `ls -m /home/genomerussia/main/fastq/unpacked/GR00/${SAMPLE}/* | sed -r "s/, /,/g" | tr -d '\n'`;
    done
