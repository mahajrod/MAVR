#!/usr/bin/env bash

PEP_DIR=$1
LABELED_PEP_DIR=$2

#SPECIES_NAME_FILE="species_names_handled.t"
mkdir -p ${LABELED_PEP_DIR}

for SPECIES_PEP_FILE in `ls ${PEP_DIR}`;
	do
	SPECIES=`basename ${SPECIES_PEP_FILE} .pep`
	echo "Labeling ${SPECIES}"
	PEP_FILE=${PEP_DIR}${SPECIES}".pep"
	LABELED_PEP_FILE=${LABELED_PEP_DIR}${SPECIES}".pep"
	sed "s/>/>${SPECIES}\./" ${PEP_FILE} > ${LABELED_PEP_FILE} 
	done




