#!/usr/bin/env bash

FUNGI_GENOMES_DIR="/home/mahajrod/Reference_genomes/fungi/ensembl/all_sm_toplevel/"
WORKDIR="/home/mahajrod/Genetics/Projects/nxf/fungi/exonerate/"
MEX67_PEP_FILE="/home/mahajrod/Genetics/Projects/nxf/fungi/MEX67_Scer.fasta"

for GENOME in `ls ${FUNGI_GENOMES_DIR}`;
    do
    GENOME_PATH=${FUNGI_GENOMES_DIR}${GENOME}
    EXONERATE_OUTPUT=${WORKDIR}${GENOME}.MEX67.exonerate.out
    exonerate  --model protein2genome --showalignment --showquerygff --showtargetgff \
               -n 100 -q ${MEX67_PEP_FILE} -t ${GENOME_PATH} > ${EXONERATE_OUTPUT}
    donehttps://scholar.google.ru/citations?hl=en&user=Ri-YMYgAAAAJ&view_op=list_works&sortby=pubdate
