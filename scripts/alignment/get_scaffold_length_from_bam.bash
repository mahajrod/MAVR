#!/usr/bin/env bash

samtools view -H $1 | grep -P "^@SQ" | awk '{printf "%s\t%s\n", $2, $3}' | sed 's/SN://;s/LN://' | sort -k2,2nr
