#!/usr/bin/env bash
WORKDIR="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/polymorphisms/"
ALIGNMENTS_DIR=${WORKDIR}"alignments/"

SUBSTITUTION_FILE=${WORKDIR}"LAN210_substitutions.t"

EXTRACTION_SCRIPT=~/Genetics/MAVR/scripts/multiple_alignment/extract_substitutions.py

cd ${WORKDIR}

echo -e ".\tLAN210_v0.10m" > ${SUBSTITUTION_FILE}

for ALIGNMENT_FILE in `ls ${ALIGNMENTS_DIR}/*`;
    do
    GENE=`basename ${ALIGNMENT_FILE} .fasta`
    #echo -en "${GENE}\t" >> ${SUBSTITUTION_FILE}
    ${EXTRACTION_SCRIPT} -i ${ALIGNMENT_FILE} -r S288C_${GENE} | sed 's/.*_//' >> ${SUBSTITUTION_FILE}
    #echo -n "\n" >> ${SUBSTITUTION_FILE}
    done
