#!/usr/bin/env bash

MAVR_DIR=$1
DATA_PREFIX=$2
PFAM_HMMSCAN_DOM_TBLOUT_FILE=$3
SWISSPROT_BLAST_TAB_HIT_FILE=$4

STAT_FILE=${DATA_PREFIX}.stat
COMPLETE_PEP_DIR="complete/"
COMPLETE_PEP_100_DIR="complete.100+/"
PEP_100_DIR="100+/"
mkdir -p ${COMPLETE_PEP_DIR} ${COMPLETE_PEP_100_DIR} ${PEP_100_DIR}

${MAVR_DIR}/scripts/sequence/histogram_length.py -w 50 -e png,svg \
                                                 -i ${DATA_PREFIX}.pep \
                                                 -o ${DATA_PREFIX}.pep

${MAVR_DIR}/scripts/transcriptome/annotation/extract_complete_proteins_from_transdecoder_output.py -i ${DATA_PREFIX}.pep \
                                                                                                   -o ${DATA_PREFIX}.complete.pep \
                                                                                                   -d ${DATA_PREFIX}.complete.pep.ids

${MAVR_DIR}/scripts/sequence/extract_sequences_by_length.py -n 100 \
                                                            -i ${DATA_PREFIX}.complete.pep \
                                                            -o ${DATA_PREFIX}.complete.100+.pep \
                                                            -d ${DATA_PREFIX}.complete.100+.pep.ids

${MAVR_DIR}/scripts/sequence/extract_sequences_by_length.py -n 100 \
                                                            -i ${DATA_PREFIX}.pep \
                                                            -o ${DATA_PREFIX}.100+.pep \
                                                            -d ${DATA_PREFIX}.100+.pep.ids
for DATA in ${DATA_PREFIX}.complete ${DATA_PREFIX}.100+ ${DATA_PREFIX}.complete.100+;
    do
    ${MAVR_DIR}/scripts/sequence/histogram_length.py -w 50 -e png,svg \
                                                     -i ${DATA}.pep \
                                                     -o ${DATA}.pep
    ${MAVR_DIR}/scripts/hmmer3/extract_hits_by_query_ids.py -i ${PFAM_HMMSCAN_DOM_TBLOUT_FILE} -f domtblout \
                                                            -d ${DATA}.pep.ids \
                                                            -o ${DATA}.pfam.domtblout
    ${MAVR_DIR}/scripts/blast/filter_blast_output_by_ids.py -m query -f blast-tab \
                                                            -w ${DATA}.pep.ids  \
                                                            -i ${SWISSPROT_BLAST_TAB_HIT_FILE} \
                                                            -o ${DATA}.swissprot.hits
    awk -F'|' '{printf "%s|%s\n", $1, $2}' ${DATA}.pep.ids > ${DATA}.mRNA.ids
    ${MAVR_DIR}/scripts/sequence/extract_sequences_by_ids.py -i ${DATA_PREFIX}.cds \
                                                             -o ${DATA}.cds \
                                                             -d ${DATA}.pep.ids
    done

#${MAVR_DIR}/scripts/sequence/histogram_length.py -w 50 -e png,svg \
#                                                 -i ${DATA_PREFIX}.complete.pep \
#                                                 -o ${DATA_PREFIX}.complete.pep
#
#${MAVR_DIR}/scripts/sequence/histogram_length.py -w 50 -e png,svg \
#                                                 -i ${DATA_PREFIX}.complete.100+.pep \
#                                                 -o ${DATA_PREFIX}.complete.100+.pep
#${MAVR_DIR}/scripts/sequence/histogram_length.py -w 50 -e png,svg \
#                                                 -i ${DATA_PREFIX}.100+.pep \
#                                                 -o ${DATA_PREFIX}.100+.pep

#${MAVR_DIR}/scripts/hmmer3/extract_hits_by_query_ids.py -i ${PFAM_HMMSCAN_DOM_TBLOUT_FILE} -f domtblout \
#                                                        -d ${DATA_PREFIX}.100+.pep.ids \
#                                                        -o ${DATA_PREFIX}.100+.pfam.domtblout
#${MAVR_DIR}/scripts/hmmer3/extract_hits_by_query_ids.py -i ${PFAM_HMMSCAN_DOM_TBLOUT_FILE} -f domtblout \
#                                                        -d ${DATA_PREFIX}.complete.pep.ids \
#                                                        -o ${DATA_PREFIX}.complete.pfam.domtblout
#${MAVR_DIR}/scripts/hmmer3/extract_hits_by_query_ids.py -i ${PFAM_HMMSCAN_DOM_TBLOUT_FILE} -f domtblout \
#                                                        -d ${DATA_PREFIX}.complete.100+.pep.ids \
#                                                        -o ${DATA_PREFIX}.complete.100+.pfam.domtblout

#${MAVR_DIR}/scripts/blast/filter_blast_output_by_ids.py -m query -f blast-tab \
#                                                        -w ${DATA_PREFIX}.100+.pep.ids  \
#                                                        -i ${SWISSPROT_BLAST_TAB_HIT_FILE} \
#                                                        -o ${DATA_PREFIX}.100+.swissprot.hits
#${MAVR_DIR}/scripts/blast/filter_blast_output_by_ids.py -m query -f blast-tab \
#                                                        -w ${DATA_PREFIX}.complete.pep.ids  \
#                                                        -i ${SWISSPROT_BLAST_TAB_HIT_FILE} \
#                                                        -o ${DATA_PREFIX}.complete.swissprot.hits
#${MAVR_DIR}/scripts/blast/filter_blast_output_by_ids.py -m query -f blast-tab \
#                                                        -w ${DATA_PREFIX}.complete.100+.pep.ids  \
#                                                        -i ${SWISSPROT_BLAST_TAB_HIT_FILE} \
#                                                        -o ${DATA_PREFIX}.complete.100+.swissprot.hits

awk -F'|' '{printf "%s|%s\n", $1, $2}' ${DATA_PREFIX}.complete.100+.pep.ids > ${DATA_PREFIX}.complete.100+.mRNA.ids
awk -F'|' '{printf "%s|%s\n", $1, $2}' ${DATA_PREFIX}.complete.pep.ids > ${DATA_PREFIX}.complete.mRNA.ids
awk -F'|' '{printf "%s|%s\n", $1, $2}' ${DATA_PREFIX}.100+.pep.ids > ${DATA_PREFIX}.100+.mRNA.ids


ALL_PEP_NUMBER=`grep -cP "^>" ${DATA_PREFIX}.pep`
ALL_PEP_LONGER_100_NUMBER=`grep -cP "^>" ${DATA_PREFIX}.100+.pep`
COMPLETE_PEP_NUMBER=`grep -cP "^>"  ${DATA_PREFIX}.complete.pep`
COMPLETE_PEP_LONGER_100_NUMBER=`grep -cP "^>" ${DATA_PREFIX}.complete.100+.pep`

echo -e "Totally peptides\t${ALL_PEP_NUMBER}" > ${STAT_FILE}
echo -e "Peptides(100+ AA)\t${ALL_PEP_LONGER_100_NUMBER}" >> ${STAT_FILE}
echo -e "Complete peptides\t${COMPLETE_PEP_NUMBER}" >> ${STAT_FILE}
echo -e "Complete peptides(100+ AA)\t${COMPLETE_PEP_LONGER_100_NUMBER}" >> ${STAT_FILE}


#${MAVR_DIR}/scripts/sequence/extract_sequences_by_ids.py -i ${DATA_PREFIX}.cds \
#                                                         -o ${DATA_PREFIX}.complete.cds \
#                                                         -d ${DATA_PREFIX}.complete.pep.ids
#
#${MAVR_DIR}/scripts/sequence/extract_sequences_by_ids.py -i ${DATA_PREFIX}.cds \
#                                                         -o ${DATA_PREFIX}.100+.cds \
#                                                         -d ${DATA_PREFIX}.100+.pep.ids
#
#${MAVR_DIR}/scripts/sequence/extract_sequences_by_ids.py -i ${DATA_PREFIX}.cds \
#                                                         -o ${DATA_PREFIX}.complete.100+.cds \
#                                                         -d ${DATA_PREFIX}.complete.100+.pep.ids

#grep -f ${DATA_PREFIX}.complete.100+.mRNA.ids ${DATA_PREFIX}.gff3 > ${DATA_PREFIX}.complete.100+.gff3
#grep -f ${DATA_PREFIX}.complete.mRNA.ids ${DATA_PREFIX}.gff3 > ${DATA_PREFIX}.complete.gff3
#grep -f ${DATA_PREFIX}.100+.mRNA.ids ${DATA_PREFIX}.gff3 > ${DATA_PREFIX}.100+.gff3

mv -f *.complete.100+.* ${COMPLETE_PEP_100_DIR}
mv -f *.complete.* ${COMPLETE_PEP_DIR}
mv -f *.100+.* ${PEP_100_DIR}