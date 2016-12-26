#!/usr/bin/env bash

#$1 - directory with assemblies
#$2 - output_file_with_statistics

wc -l $1 | sort -k1nr | awk 'BEGIN {printf "Genome_id\tNumber_of_assemblies\n"};NR>1 {split($2,ID_ARRAY,".");printf "%i\t%s\n", ID_ARRAY[1], $1-1}' > $2
