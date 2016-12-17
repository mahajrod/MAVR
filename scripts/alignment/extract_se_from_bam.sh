#!/usr/bin/env bash

samtools view $1 | awk -F'\t' 'BEGIN {A=1}; {printf "@READ%i_%s\n%s\n+\n%s\n",A, $1,$10,$11; A+=1}' > $2
