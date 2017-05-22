#!/usr/bin/env bash

for SP in m54099_170411_113141 m54099_170411_173102 m54099_170412_113507; do ~/soft/MAVR/scripts/alignment/get_read_names_from_sam.py -i ${SP}.bam -o ${SP}.read_names.ids; sed 's/\/ccs$//' ${SP}.read_names.ids > ${SP}.read_names.prefix.ids;  samtools view -h ${SP}.subreads.bam | ~/soft/MAVR/scripts/alignment/get_reads_by_name.py -r ${SP}.read_names.prefix.ids -m include -c partial | samtools view -hb -o - ${SP}.subreads_for_ccs.bam &  samtools view -h ${SP}.subreads.bam | ~/soft/MAVR/scripts/alignment/get_reads_by_name.py -r ${SP}.read_names.prefix.ids -m remove -c partial | samtools view -hb -o - ${SP}.subreads.no_ccs.bam - & done

for SP in m54099_170411_113141 m54099_170411_173102 m54099_170412_113507;
do ~/soft/MAVR/scripts/alignment/get_read_names_from_sam.py -i ${SP}.bam -o ${SP}.read_names.ids;
sed 's/\/ccs$//' ${SP}.read_names.ids > ${SP}.read_names.prefix.ids;
samtools view -h ${SP}.subreads.bam | ~/soft/MAVR/scripts/alignment/get_reads_by_name.py -r ${SP}.read_names.prefix.ids -m include -c partial | samtools view -b -o ${SP}.subreads_for_ccs.bam - &
samtools view -h ${SP}.subreads.bam | ~/soft/MAVR/scripts/alignment/get_reads_by_name.py -r ${SP}.read_names.prefix.ids -m remove -c partial | samtools view -b -o ${SP}.subreads.no_ccs.bam - &
done