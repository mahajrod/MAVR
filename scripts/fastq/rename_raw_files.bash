#!/usr/bin/env bash

for SP in $0; do mkdir ${SP}; mv ${SP}_* ${SP}; done

for SP in $0; do mv ${SP}/${SP}_*R1*fastq.gz ${SP}/${SP}_1.fastq.gz; done
for SP in $0; do mv ${SP}/${SP}_*R2*fastq.gz ${SP}/${SP}_2.fastq.gz; done
