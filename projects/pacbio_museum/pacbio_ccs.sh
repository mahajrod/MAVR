#!/usr/bin/env bash

for SP in m54099_170411_113141 m54099_170411_173102 m54099_170412_113507; do  /nfs/mfnstore-lin/export/mfn-genom-1/KLIVER/smrtlink/soft/smrtcmds/bin/ccs --maxLength=15000  --minLength=100 --numThreads=64 --minPasses=1 --reportFile=${SP}.report ${SP}.subreads.bam ${SP}.bam; done

for SP in m54099_170411_113141 m54099_170411_173102 m54099_170412_113507; do ~/soft/MAVR/scripts/alignment/get_read_names_from_sam.py -i ${SP}.bam -o ${SP}.read_names.ids; sed 's/\/ccs$//' ${SP}.read_names.ids > ${SP}.read_names.prefix.ids;  samtools view -h ${SP}.subreads.bam | ~/soft/MAVR/scripts/alignment/get_reads_by_name.py -r ${SP}.read_names.prefix.ids -m include -c partial | samtools view -hb -o - ${SP}.subreads_for_ccs.bam &  samtools view -h ${SP}.subreads.bam | ~/soft/MAVR/scripts/alignment/get_reads_by_name.py -r ${SP}.read_names.prefix.ids -m remove -c partial | samtools view -hb -o - ${SP}.subreads.no_ccs.bam - & done

for SP in m54099_170411_113141 m54099_170411_173102 m54099_170412_113507;
    do ~/soft/MAVR/scripts/alignment/get_read_names_from_sam.py -i ${SP}.bam -o ${SP}.read_names.ids;
    sed 's/\/ccs$//' ${SP}.read_names.ids > ${SP}.read_names.prefix.ids;
    samtools view -h ${SP}.subreads.bam | ~/soft/MAVR/scripts/alignment/get_reads_by_name.py -r ${SP}.read_names.prefix.ids -m include -c partial | samtools view -b -o ${SP}.subreads_for_ccs.bam - &
    samtools view -h ${SP}.subreads.bam | ~/soft/MAVR/scripts/alignment/get_reads_by_name.py -r ${SP}.read_names.prefix.ids -m remove -c partial | samtools view -b -o ${SP}.subreads.no_ccs.bam - &
    done


#ccs call
for SP in m54099_170515_135132 m54099_170515_195104 m54099_170517_145238 m54099_170517_205155 m54099_170522_131017 m54099_170522_190936 m54099_170523_100047 m54099_170523_160014 m54099_170523_221222; do  /nfs/mfnstore-lin/export/mfn-genom-1/KLIVER/smrtlink/soft/smrtcmds/bin/ccs --maxLength=15000  --minLength=100 --numThreads=64 --minPasses=1 --reportFile=${SP}.report ${SP}.subreads.bam ${SP}.bam; done

#split collapsed and not collapsed reads
for SP in m54099_170515_135132 m54099_170515_195104 m54099_170517_145238 m54099_170517_205155 m54099_170522_131017 m54099_170522_190936 m54099_170523_100047 m54099_170523_160014 m54099_170523_221222;
    do ~/soft/MAVR/scripts/alignment/get_read_names_from_sam.py -i ${SP}.bam -o ${SP}.read_names.ids;
    sed 's/\/ccs$//' ${SP}.read_names.ids > ${SP}.read_names.prefix.ids;
    samtools view -h ${SP}.subreads.bam | ~/soft/MAVR/scripts/alignment/get_reads_by_name.py -r ${SP}.read_names.prefix.ids -m include -c partial | samtools view -b -o ${SP}.subreads_for_ccs.bam - &
    samtools view -h ${SP}.subreads.bam | ~/soft/MAVR/scripts/alignment/get_reads_by_name.py -r ${SP}.read_names.prefix.ids -m remove -c partial | samtools view -b -o ${SP}.subreads.no_ccs.bam - &
    done