#!/usr/bin/env bash

cd /mnt/peru/skliver/cat_comparison/augustus

for SP in *;
    do
    echo ${SP};
    ~/Soft/MAVR/scripts/sequence_clusters/extract_clusters_by_element_ids.py -i ${SP}/${SP}.ortho.orthologs \
                                                                             -d ${SP}/supported_by_db_and_hints/${SP}.augustus.supported.transcripts.swissprot_or_pfam_and_hints.longest_pep.evidence.ids \
                                                                             -o ${SP}/supported_by_db_and_hints/${SP}.augustus.supported.transcripts.swissprot_or_pfam_and_hints.longest_pep.evidence.orthologs;
    grep -f /mnt/guatemala/skliver/data/EggOG/maNOG/maNOG.olfactory.ids ${SP}/supported_by_db_and_hints/${SP}.augustus.supported.transcripts.swissprot_or_pfam_and_hints.longest_pep.evidence.orthologs > ${SP}/supported_by_db_and_hints/${SP}.augustus.supported.transcripts.swissprot_or_pfam_and_hints.longest_pep.evidence.olfactory.fam;
    ~/Soft/MAVR/scripts/sequence_clusters/treefam/extract_proteins_from_selected_families.py -i /mnt/guatemala/skliver/data/EggOG/maNOG/maNOG.olfactory.ids \
                                                                                             -f ${SP}/supported_by_db_and_hints/${SP}.augustus.supported.transcripts.swissprot_or_pfam_and_hints.longest_pep.evidence.olfactory.fam \
                                                                                             -p ${SP}/supported_by_db_and_hints/${SP}.augustus.supported.transcripts.swissprot_or_pfam_and_hints.longest_pep.pep \
                                                                                             -d ${SP}/supported_by_db_and_hints/olfactory;
    done
