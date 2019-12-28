#!/usr/bin/env bash
#TODO: adjust script. In fact it is a template
for SP in $1; do mkdir -p ${SP}; cat ${SP}*ids | xargs -P $2 -I SRA fasterq-dump --progress --mem 10000MB -O ${SP} SRA & done
