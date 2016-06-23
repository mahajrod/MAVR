#!/usr/bin/env bash

FAM_DIR=$1
LABELED_FAM_DIR=$2

#SPECIES_NAME_FILE="species_names_handled.t"
mkdir -p ${LABELED_FAM_DIR}

for SPECIES_FAM_FILE in `ls ${FAM_DIR}`;
	do
	SPECIES=`basename ${SPECIES_FAM_FILE} .orthologs`
	echo "Labeling ${SPECIES}"
	FAM_FILE=${FAM_DIR}/${SPECIES}".orthologs"
	LABELED_FAM_FILE=${LABELED_FAM_DIR}/${SPECIES}".orthologs"
	sed "s/\t/\t${SPECIES}\@/;s/,/,${SPECIES}\@/g" ${FAM_FILE} > ${LABELED_FAM_FILE}
	done




