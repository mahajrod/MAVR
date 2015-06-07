#!/usr/bin/env bash

PEP_DIR="pep/"
LABELED_PEP_DIR="labeled_pep/"

SPECIES_NAME_FILE="species_names_handled.t"
mkdir -p ${LABELED_PEP_DIR}

for SPECIES in `cat ${SPECIES_NAME_FILE}`;
	do
	echo "Labeling ${SPECIES}"
	PEP_FILE=${PEP_DIR}${SPECIES}".pep"
	LABELED_PEP_FILE=${LABELED_PEP_DIR}${SPECIES}".pep"
	sed "s/>/>${SPECIES}\./" ${PEP_FILE} > ${LABELED_PEP_FILE} 
	done




