#!/usr/bin/env bash

#skliver@supermicro:~/workdir/single_gene_families/species_set_1$

~/soft/MAVR/scripts/treefam/label_pep_with_species_name.sh pep/ labeled_pep/
~/soft/MAVR/scripts/treefam/label_fam_file_with_species_name.sh fam/ labeled_fam

for SPECIES in acinonyx_jubatus ailuropoda_melanoleuca bos_taurus canis_familiaris equus_caballus felis_catus homo_sapiens loxodonta_africana macaca_mulatta monodelphis_domestica mus_musculus mustela_putorius_furo myotis_lucifugus odobenus_rosmarus oryctolagus_cuniculus panthera_tigris pteropus_alecto pusa_sibirica rattus_norvegicus sarcophilus_harrisii sus_scrofa ursus_maritimus;
    do
    ~/soft/MAVR/scripts/treefam/extract_proteins_from_selected_families.py -i single_gene_families_species_set1.ids \
                                                                           -f labeled_fam/${SPECIES}.fam \
                                                                           -p labeled_pep/${SPECIES}.pep \
                                                                           -c -d fam_pep/ \
                                                                           -o ${SPECIES};
    done

mkdir combined_fam_pep

for FAMILY in `ls ./fam_pep/ `;
    do
    echo $FAMILY;
    cat ./fam_pep/${FAMILY}/* > combined_fam_pep/${FAMILY}.pep;
    done

~/soft/MAVR/scripts/multiple_alignment/parallel_mafft.py -p 16 -i combined_fam_pep/ -o fam_alignments/