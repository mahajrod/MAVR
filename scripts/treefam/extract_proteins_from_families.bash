#!/usr/bin/env bash

THREADS=4

CURRENT_DIR=`pwd`
WORK_DIR="family_proteins/"
PROT_IDS_SUBDIR="prot_ids/"
GENE_IDS_SUBDIR="gene_ids/"
PROT_FAM_SUBDIR="prot_fam/"
GENE_FAM_SUBDIR="gene_fam/"
EMF_SUBDIR="emf/"

PROT_SEQ_SUBDIR="pep/"
CDS_SEQ_SUBDIR="cds/"
ALL_PROTEIN_DIR=${CURRENT_DIR}"/treefam9_species/pep/"
ALL_CDS_DIR=${CURRENT_DIR}"/treefam9_species/cds/"
TREEFAM_FAM_DATA_DIR=${CURRENT_DIR}"/treefam_family_data/"

SPECIES_FILE=${CURRENT_DIR}"/species_names.t"

mkdir -p ${WORK_DIR}
cd ${WORK_DIR}

mkdir -p ${PROT_IDS_SUBDIR} ${GENE_IDS_SUBDIR} ${PROT_FAM_SUBDIR} ${GENE_FAM_SUBDIR} ${EMF_SUBDIR} ${PROT_SEQ_SUBDIR} ${CDS_SEQ_SUBDIR}

for SPECIES in `cat ${SPECIES_FILE}`;
	do
	
	echo "Handling ${SPECIES}..."
	
	ALL_PROT_FILE=${ALL_PROTEIN_DIR}${SPECIES}".fa"
	ALL_CDS_FILE=${ALL_CDS_DIR}${SPECIES}".cds.fa"
	PROT_IDS_FILE=${SPECIES}"_in_treefam_families_protein.ids"
	GENE_IDS_FILE=${SPECIES}"_in_treefam_families_genes.ids"
	PROT_FAM_FILE=${SPECIES}"_proteins.fam"
	GENE_FAM_FILE=${SPECIES}"_genes.fam"
	PROT_SEQ_FILE=${PROT_SEQ_SUBDIR}${SPECIES}".pep" 
	CDS_SEQ_FILE=${CDS_SEQ_SUBDIR}${SPECIES}".cds"

	~/Genetics/MAVR/scripts/treefam/extract_families_of_species.py -p ${THREADS} -s ${SPECIES} -i ${TREEFAM_FAM_DATA_DIR}
	~/Genetics/MAVR/scripts/sequence/extract_sequences_by_ids.py -i ${ALL_PROT_FILE} -d ${PROT_IDS_FILE} -o ${PROT_SEQ_FILE}
	~/Genetics/MAVR/scripts/sequence/extract_sequences_by_ids.py -i ${ALL_CDS_FILE} -d ${PROT_IDS_FILE} -o ${CDS_SEQ_FILE}
	
	mv ${PROT_IDS_FILE} ${PROT_IDS_SUBDIR}
	mv ${GENE_IDS_FILE} ${GENE_IDS_SUBDIR}
	mv ${PROT_FAM_FILE} ${PROT_FAM_SUBDIR}
	mv ${GENE_FAM_FILE} ${GENE_FAM_SUBDIR}
	mv ${SPECIES} ${EMF_SUBDIR}

	done

