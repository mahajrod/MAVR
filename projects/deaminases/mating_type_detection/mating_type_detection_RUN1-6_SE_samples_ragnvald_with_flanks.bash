#!/usr/bin/env bash

source ./DESCRIPTION

BOWTIE2_INDEX=${WORKDIR}"index/MatA1_MatAlpha1_with_100_flanks"

MAVR_SCRIPTS_DIR=${MAVR_DIR}"/scripts/"
RESTORE_PAIRS_SCRIPT=${MAVR_SCRIPTS_DIR}"filter/restore_pairs.py"
MAP_SCRIPT=${MAVR_SCRIPTS_DIR}"alignment/map_reads.py"
RESTORE_PAIRS_SCRIPT=${MAVR_SCRIPTS_DIR}"filter/restore_pairs.py"

EXTRACTOR_BIN="extract"
FRAGMENTS_FILE=${WORKDIR}"kmer/MatA1_MatAlpha1_with_100_flanks_with_rev_com_33_mer.kmer"
SAMPLES_DIR="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all/"

SAMPLES=(N041-LAN210-Can-AID-Oct12-RUN4-D3 N050-LAN211-Can-HAP-NA-RUN4 N034-LAN210-Can-A1-Oct12-RUN4-D3 N039-LAN210-Can-A3G-Oct12-RUN4-D3 N036-LAN210-FOA-A1-Oct12-RUN4-D3 N051-LAN211-Can-HAP-NA-RUN4 N037-LAN210-FOA-A1-Oct12-RUN4-D3 N038-LAN210-Can-A3G-Oct12-RUN4-D3 N035-LAN210-Can-A1-Oct12-RUN4-D3 N040-LAN210-Can-AID-Oct12-RUN4-D3 N041-LAN210-Can-AID-Oct12-RUN4-D3 N050-LAN211-Can-HAP-NA-RUN4 N034-LAN210-Can-A1-Oct12-RUN4-D3 N039-LAN210-Can-A3G-Oct12-RUN4-D3 N036-LAN210-FOA-A1-Oct12-RUN4-D3 N051-LAN211-Can-HAP-NA-RUN4 N037-LAN210-FOA-A1-Oct12-RUN4-D3 N038-LAN210-Can-A3G-Oct12-RUN4-D3 N035-LAN210-Can-A1-Oct12-RUN4-D3 N040-LAN210-Can-AID-Oct12-RUN4-D3)

for SAMPLE in ${SAMPLES[@]};
    do
    cd ${WORKDIR}

    CURRENT_SAMPLE=${SAMPLE}

    mkdir -p ${CURRENT_SAMPLE}
    cd ${CURRENT_SAMPLE}
    echo ${CURRENT_SAMPLE}

    READS=(`ls ${SAMPLES_DIR}${SAMPLE}/trimmed/spades/corrected/ | grep -P ".*.fastq$"`)
    READS_FILE=${SAMPLES_DIR}${SAMPLE}/trimmed/spades/corrected/${READS}

    ${EXTRACTOR_BIN} --fragments ${FRAGMENTS_FILE} -o ./ -i ${READS_FILE}

    READ_BASENAME=`basename ${READS[0]} .fastq`

    mkdir -p all_reads

    EXTRACTED_READS_FILE=${READ_BASENAME}".filtered.fastq"
    cp ${EXTRACTED_READS_FILE} "all_reads/"${EXTRACTED_READS_FILE}

    ALL_PREFIX=./all_reads/${CURRENT_SAMPLE}"_all"


    ${MAP_SCRIPT} -i ${BOWTIE2_INDEX} -t 4 -u ${EXTRACTED_READS_FILE} \
                  -p ${ALL_PREFIX} -g -y ${ALL_PREFIX}"_coverage.bed" -z -x -f 100

    done
