#!/usr/bin/env bash
#dell
AUGUSTUS_CONFIG_PATH=/storage1/home/skliver/Soft/augustus-3.2.1/config/
GENEMARK_PATH=/storage1/home/skliver/Soft/GeneMark/gm_et_linux_64/gmes_petap/
BAMTOOLS_PATH=/storage1/home/skliver/Soft/bamtools/bin

cd /storage1/home/skliver/workdir/tiger_transcriptome/braker
~/Soft/BRAKER1/braker.pl --genome ../scaffold_merged.fa \
                         --bam ../siberian_tiger.bam \
                         --cores 15 \
                         --species=Panthera_tigris_amur
