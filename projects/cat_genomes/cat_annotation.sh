#!/usr/bin/env bash


SPECIES="felis_catus"
SPECIES_NAME="Felis catus"
PEP_DIR="/home/skliver/data/genomes/cat/v8.0/pep/"
RAW_PEP_FILE="GCF_000181335.2_Felis_catus_8.0_protein.faa"

cd ${PEP_DIR}

grep -P "^>" ${RAW_PEP_FILE} | sed 's/^>//;s/\[${SPECIES_NAME}\]//;s/ /\t/' > ${SPECIES}.pep.description

sort -t $'\t' -k 2 -k 1 ${SPECIES}.pep.description > ${SPECIES}.pep.sorted.description

#uniq ignoring first field
uniq -f 1 ${SPECIES}.pep.sorted.description > ${SPECIES}.pep.sorted.description.uniq

sed s/isoform.*// ${SPECIES}.pep.sorted.description.uniq > ${SPECIES}.pep.sorted.description.no_isoform_versions.uniq

~/soft/MAVR/scripts/collaps_synonym_strings.py -i ${SPECIES}.pep.sorted.description.no_isoform_versions.uniq \
                                               -k 1 -a 0 -o ${SPECIES}.pep.collapsed_isoforms.description



