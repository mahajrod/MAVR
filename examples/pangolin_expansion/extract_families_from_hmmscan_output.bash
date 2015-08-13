#!/usr/bin/env bash

# supermicro

MAVR_DIR=~/soft/MAVR/
SCRIPTS_DIR=${MAVR_DIR}scripts/hmmer3/
#EXTRACT_HITS_SCRIPT=${SCRIPTS_DIR}extract_hits_from_hmm_output.py
TOP_HITS_SCRIPT=${SCRIPTS_DIR}extract_top_hits_from_hmm_output.py
FAMILIES_SCRIPT=${SCRIPTS_DIR}convert_top_hits_to_families.py

#for SPECIES in manis_pentadactyla manis_javanica;
#    do
#    cd ~/workdir/pangolin/${SPECIES}/
#    $TOP_HITS_SCRIPT -i ${SPECIES}_hmmscan.report \
#                    -o ${SPECIES}_hmmscan_top_hits.tab \
#                    -f hmmer3-text -n ${SPECIES}_hmmscan_not_found.t \
#                    -g ${SPECIES}_hmmscan_not_significant.t
#
#    $FAMILIES_SCRIPT -i ${SPECIES}_hmmscan_top_hits.tab \
#                    -e -k 1 -a 0 \
#                    -o ${SPECIES}_hmmscan_families.tab
#    done

cd /home/skliver/workdir/pangolin/hmmer_expansion/

${MAVR_DIR}scripts/treefam/prepare_cafe_input_from_fam_files.py -i /home/skliver/workdir/pangolin/hmmer_expansion/prot_fam_sets \
                                                                -c 8_sp_2_pangolin_not_filtered.cafe \
                                                                -s homo_sapiens,mus_musculus,pteropus_vampyrus,bos_taurus,equus_caballus,felis_catus,ailuropoda_melanoleuca,canis_familiaris,manis_javanica,manis_pentadactyla \
                                                                -u .fam \
${MAVR_DIR}scripts/treefam/prepare_cafe_input_from_fam_files.py -i /home/skliver/workdir/pangolin/hmmer_expansion/prot_fam_sets \
                                                                -c 8_sp_2_pangolin_filtered_ambigious_fams.cafe \
                                                                -s homo_sapiens,mus_musculus,pteropus_vampyrus,bos_taurus,equus_caballus,felis_catus,ailuropoda_melanoleuca,canis_familiaris,manis_javanica,manis_pentadactyla \
                                                                -u .fam \
                                                                -w ~/data/TreeFam_new/check_ambigious/families_all_species_histo_errors_less_0.05.ids

${MAVR_DIR}scripts/treefam/prepare_cafe_input_from_fam_files.py -i /home/skliver/workdir/pangolin/hmmer_expansion/prot_fam_sets \
                                                                -c 8_sp_2_pangolin_filtered_ambigious_fams_and_fams_with_less_then_3_sp.cafe \
                                                                -s homo_sapiens,mus_musculus,pteropus_vampyrus,bos_taurus,equus_caballus,felis_catus,ailuropoda_melanoleuca,canis_familiaris,manis_javanica,manis_pentadactyla \
                                                                -u .fam \
                                                                -m 3 \
                                                                -w ~/data/TreeFam_new/check_ambigious/families_all_species_histo_errors_less_0.05.ids


