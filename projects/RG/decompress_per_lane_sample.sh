#!/usr/bin/env bash

SAMPLE="Sample_27F"; for FILE in `ls ${SAMPLE}/*.f*q.gz`; do echo "Handling ${FILE}"; mkdir -p ../fastq/${SAMPLE}; OUT_FILE=../fastq/${SAMPLE}/`basename ${FILE} .gz`; echo "Decompressing to ${OUT_FILE}"; gunzip -c ${FILE} > ${OUT_FILE}; done