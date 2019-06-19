#!/usr/bin/env bash
#TODO: adjust script. In fact it is a template
for SP in $1; do mkdir -p ${SP}; cat ${SP}*ids | xargs -P 3 -I NAME fastq-dump --split-3 -F -O ${SP} NAME & done
