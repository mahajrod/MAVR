#!/usr/bin/env bash
# supermicro

WORK_DIR="/home/skliver/data/TreeFam_new/check_ambigious/"
MAVR_DIR=~/soft/MAVR/
SCRIPTS_DIR=${MAVR_DIR}scripts/
SPECIES_FAM_DIR="/home/skliver/data/TreeFam_new/check_ambigious/species_fam/"
REFERENCE_FAM_DIR="/home/skliver/data/TreeFam_new/fam/"
COMPARE_SCRIPT=${SCRIPTS_DIR}expansion/compare_clusters.py
FIND_CORRECT_FAMS_SCRIPT=${SCRIPTS_DIR}expansion_hmm/find_correctly_assembled_families_in_all_species.py
GATHER_STAT_SCRIPT=${SCRIPTS_DIR}expansion/gather_statistics_from_comparison.py
ERRORS_DIST_SCRIPT=${SCRIPTS_DIR}expansion_hmm/errors_distribution.py

SPECIES=(rattus_norvegicus macaca_mulatta echinops_telfairi ornithorhynchus_anatinus monodelphis_domestica tupaia_belangeri erinaceus_europaeus sorex_araneus microcebus_murinus pongo_abelii equus_caballus felis_catus ochotona_princeps cavia_porcellus choloepus_hoffmanni procavia_capensis tursiops_truncatus tarsius_syrichta dipodomys_ordii vicugna_pacos pteropus_vampyrus dasypus_novemcinctus homo_sapiens macropus_eugenii loxodonta_africana oryctolagus_cuniculus ailuropoda_melanoleuca nomascus_leucogenys callithrix_jacchus myotis_lucifugus sarcophilus_harrisii bos_taurus gorilla_gorilla otolemur_garnettii pan_troglodytes ictidomys_tridecemlineatus sus_scrofa mus_musculus canis_familiaris mustela_putorius_furo)

cd ${WORK_DIR}
mkdir -p compare_dir
for SP in ${SPECIES[@]};
    do

    mkdir -p compare_dir/${SP}
    echo "comparing families of ${SP}"
    ${COMPARE_SCRIPT} -r ${REFERENCE_FAM_DIR}${SP}.fam -c ${SPECIES_FAM_DIR}${SP}.fam -n -o compare_dir/${SP}

    done

SPECIES_STR=`echo ${SPECIES[@]} | tr " " ","`

${GATHER_STAT_SCRIPT} -s ${SPECIES_STR} -d compare_dir/ -o mammalia_hmmscan_families.stat
${FIND_CORRECT_FAMS_SCRIPT} -s ${SPECIES_STR} -d compare_dir/
${ERRORS_DIST_SCRIPT} -e -i families_all_species.t | awk 'NR == 1; NR > 1 {print $0 | "sort -nr -k 5"}' > families_all_species_histo.stat

awk 'NR == 1; NR > 1 {print $0 | if ($5 <= 0.05) {print $0} }' families_all_species_histo.stat
awk 'NR == 1; NR > 1 {if ($5 <= 0.05) {print $0}; }' families_all_species_histo.stat > families_all_species_histo_errors_less_0.05.stat

#this families may be used in following expansion analysis
awk ' NR > 1 {if ($5 <= 0.05) {print $1}; }' families_all_species_histo.stat > families_all_species_histo_errors_less_0.05.ids