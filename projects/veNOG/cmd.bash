#!/usr/bin/env bash

~/Soft/MAVR/scripts/sequence_clusters/eggnog/extract_proteins_from_alignments.py -d  veNOG_raw_algs/ -o veNOG_proteins

~/Soft/MAVR/scripts/sequence_clusters/eggnog/split_proteins_per_species.py -i  veNOG_proteins -o veNOG_proteins_per_specie

cd /mnt/guatemala/skliver/data/EggOG/veNOG/veNOG_proteins_per_species
ls | sed s/.pep// > ../veNOG.species.taxa.id

~/Soft/MAVR/scripts/sequence_clusters/eggnog/rename_files_from_taxaid_to_species_name.py -i veNOG_proteins_per_species/ -d veNOG.species.replaced_spaces.species -o veNOG_proteins_species

~/Soft/MAVR/scripts/entrez/get_taxonomy.py -i veNOG.species.taxa.ids  -o veNOG.species.taxonomy -a mahajrod@gmail.com -t id

awk -F'\t' '{printf "%s\t%s\n", $1, $3}' veNOG.species.taxonomy > veNOG.species.ids.syn
sed 's/\ /_/' veNOG.species.ids.syn > veNOG.species.ids.no_space.syn



~/Soft/MAVR/scripts/sequence_clusters/eggnog/convert_members_tsv_to_fam.py -i veNOG.members.tsv -o veNOG.members.fam

~/Soft/MAVR/scripts/sequence_clusters/split_clusters_by_element_labels.py -i veNOG.members.fam -s 10020,10090,10116,10141,132908,13616,13735,28377,30608,30611,31033,37347,43179,59463,59729,61853,69293,7757,7897,7955,8049,8083,8090,8128,8364,9031,9103,9258,9305,9315,9361,9371,9478,9483,9544,9593,9598,9601,9606,9615,9646,9669,9685,9739,9785,9796,9813,9823,9913,9986,99883 -d fam/ -e .

~/Soft/MAVR/scripts/sequence_clusters/eggnog/rename_files_from_taxaid_to_species_name.py -i fam/ -d veNOG.species.replaced_spaces.species -o fam_latin/


cd /mnt/guatemala/skliver/white_shark_project/selection
mkdir labeled_fam~/Soft/MAVR/scripts/sequence_clusters/expansion/prepare_cafe_input.py -i 9sp.fam  -c 9sp.cafe -s astyanax_mexicanus,callorhinchus_milii,carcharodon_carcharias,danio_rerio,latimeria_chalumnae,lepisosteus_oculatus,oreochromis_niloticus,poecilia_formosa,rhincodon_typus -e @

for FILE in `ls labeled_fam`; do sed 's/veNOG.//;s/.meta_raw//' fam/${FILE} > corrected_fam/${FILE}; done
for FILE in `ls fam/`; do ~/Soft/MAVR/scripts/sequence_clusters/label_cluster_elements.py -i corrected_fam/${FILE} -o labeled_fam/${FILE}; done



~/Soft/MAVR/scripts/sequence_clusters/merge_fam_files.py -i labeled_fam/ -o 9sp.fam

~/Soft/MAVR/scripts/sequence_clusters/expansion/prepare_cafe_input.py -i 9sp.fam  -c 9sp.cafe -s astyanax_mexicanus,callorhinchus_milii,carcharodon_carcharias,danio_rerio,latimeria_chalumnae,lepisosteus_oculatus,oreochromis_niloticus,poecilia_formosa,rhincodon_typus -e @

~/Soft/MAVR/scripts/sequence_clusters/extract_single_copy_clusters.py -i labeled_fam/ -o 9sp.monoclusters.fam

# Was found 2011 single-copy clusters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

mkdir labeled_pep

for FILE in `ls pep/`; do ~/Soft/MAVR/scripts/sequence/label_sequences.py -i pep/${FILE} -l `basename ${FILE} .pep` -o labeled_pep/${FILE}; done

~/.local/bin/extract_sequences_from_selected_clusters.py -f 9sp.monoclusters.fam  -p labeled_pep/ -d monoclusters

~/Soft/MAVR/scripts/multiple_alignment/remove_columns_with_gaps.py -i monoclusters_alignment -n 5 -o monoclusters_alignment_max_5_gaps/  -v

for FILE in `ls monoclusters_alignment_max_5_gaps`; do ~/Soft/MAVR/scripts/sequence/leave_only_label.sh monoclusters_alignment_max_5_gaps/${FILE}  monoclusters_alignment_max_5_gaps_only_labels/${FILE}; done

~/Soft/MAVR/scripts/multiple_alignment/merge_alignments.py -i monoclusters_alignment_max_5_gaps_only_labels/ -o 9sp_monoclusters_max5_gaps.alignment.fasta -c 9sp_monoclusters_max5_gaps.alignment.coords

 ~/Soft/raxml/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3 -f a -T 26 -# 1000 -n 9sp_shark -m PROTGAMMAAUTO -s 9sp_monoclusters_max5_gaps.alignment.no_custom_aa.fasta -x 12345 -p 1234567