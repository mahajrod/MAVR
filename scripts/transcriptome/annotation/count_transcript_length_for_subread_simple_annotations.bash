#!/usr/bin/env bash


#awk -F'\t' 'NR == 1 {LEN=$2; PREV=$1}; NR > 1{if ($1 != PREV){printf "%s\t%i\n", PREV, LEN; PREV=$1;LEN=$2} else {LEN+=$2}} END {printf "%s\t%i\n", PREV, LEN}' hg38_RefSeq_exon.ncbi_names.exon_len > hg38_RefSeq_exon.ncbi_names.transcript_len


awk -F'\t' 'NR >1 {printf "%s\t%i\n", $1, $4-$3+1}' $1 | awk -F'\t' 'NR == 1 {LEN=$2; PREV=$1}; NR > 1{if ($1 != PREV){printf "%s\t%i\n", PREV, LEN; PREV=$1;LEN=$2} else {LEN+=$2}} END {printf "%s\t%i\n", PREV, LEN}'
