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

SPECIES=(rattus_norvegicus macaca_mulatta echinops_telfairi ornithorhynchus_anatinus monodelphis_domestica tupaia_belangeri erinaceus_europaeus sorex_araneus microcebus_murinus pongo_abelii equus_caballus felis_catus ochotona_princeps cavia_porcellus choloepus_hoffmanni procavia_capensis tursiops_truncatus tarsius_syrichta dipodomys_ordii vicugna_pacos pteropus_vampyrus dasypus_novemcinctus homo_sapiens macropus_eugenii loxodonta_africana oryctolagus_cuniculus ailuropoda_melanoleuca nomascus_leucogenys callithrix_jacchus myotis_lucifugus sarcophilus_harrisii bos_taurus gorilla_gorilla otolemur_garnettii pan_troglodytes ictidomys_tridecemlineatus sus_scrofa mus_musculus canis_familiaris mustela_putorius_furo)

cd ${WORK_DIR}
mkdir -p compare_dir
for species in ${SPECIES[@]};
    do

    mkdir -p compare_dir/${SPECIES}
    echo "comparing families of ${species}"
    ${COMPARE_SCRIPT} -r ${REFERENCE_FAM_DIR}${SPECIES}.fam -c ${SPECIES_FAM_DIR}${SPECIES}.fam -n -o compare_dir/${SPECIES}

    done

SPECIES_STR=`echo ${SPECIES[@]} | tr " " ","`
${GATHER_STAT_SCRIPT} -s ${SCRIPTS_DIR} -d compare_dir/ -o mammalia_hmmscan_families.stat
${FIND_CORRECT_FAMS_SCRIPT} -s ${SPECIES_STR} -d compare_dir/