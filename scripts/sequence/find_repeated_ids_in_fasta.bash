#!/usr/bin/env bash

TOTAL_IDS=`grep -cP "^>" $1`
UNIQUE_IDS=`grep -P "^>" $1 | sed "s/ .*//" | sort | uniq | wc -l`

echo "Totaly ids: ${TOTAL_IDS}"
echo "Unique ids: ${UNIQUE_IDS}"

echo "Repeated ids:"
grep -P "^>" $1 | sed "s/ .*//" | sort | uniq -c  | awk '{if ($1 > 1) print $2}'
