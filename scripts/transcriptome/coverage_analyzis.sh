#!/usr/bin/env bash
SAMPLE="merged.alignment.sorted.q20"

ALIGNMENT_FILE="${SAMPLE}.bam"
COVERAGE_FILE="${SAMPLE}.coverage"
COLLAPSED_COVERAGE_FILE="${SAMPLE}.collapsed.coverage"
STAT_FILE="${SAMPLE}.collapsed.coverage.stats"

COVERAGE100_ID_FILE="${SAMPLE}.coverage100.ids"
COVERAGE90_ID_FILE="${SAMPLE}.coverage90.ids"
COVERAGE50_ID_FILE="${SAMPLE}.coverage50.ids"

bedtools genomecov -ibam ${FILTERED_ALIGNMENT_FILE} -d > ${COVERAGE_FILE}

awk -F'\t' 'BEGIN {SCAF=""; LEN=""; COV=""} {if (($1 != SCAF)) {if (NR > 1) {printf "%s\t%s\t%s\n",SCAF,LEN, COV}; SCAF=$1; LEN=$2; COV=$3} else {LEN=$2; COV=COV","$3}} ; END {printf "%s\t%s\t%s\n",SCAF,LEN, COV}'  ${COVERAGE_FILE} > ${COLLAPSED_COVERAGE_FILE}

~/soft/MAVR/scripts/alignment/analize_collapsed_coverage_file.py -i ${COLLAPSED_COVERAGE_FILE} -o ${STAT_FILE}

COVERAGE100=`awk -F'\t' 'NR>1 {if($9==0) print $1}' ${STAT_FILE} | tee ${COVERAGE100_ID_FILE} | wc -l`

COVERAGE90=`awk -F'\t' 'NR>1 {if($7/$2 >= 0.9) print $1}' ${STAT_FILE} | tee ${COVERAGE90_ID_FILE} | wc -l`

COVERAGE50=`awk -F'\t' 'NR>1 {if($7/$2 >= 0.5) print $1}' ${STAT_FILE} | tee ${COVERAGE50_ID_FILE} | wc -l`
