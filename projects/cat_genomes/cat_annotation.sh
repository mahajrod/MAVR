#!/usr/bin/env bash


cd /home/skliver/data/genomes/cat/v8.0/pep/

grep -P "^>" GCF_000181335.2_Felis_catus_8.0_protein.faa | sed 's/^>//;s/\[Felis catus\]//;s/ /\t/' > felis_catus.pep.description

sort -t $'\t' -k 2 -k 1 felis_catus.pep.description > felis_catus.pep.sorted.description

#uniq ignoring first field
uniq -f 1 felis_catus.pep.sorted.description > felis_catus.pep.sorted.description.uniq

sed s/isoform.*// felis_catus.pep.sorted.description.uniq > felis_catus.pep.sorted.description.no_isoform_versions.uniq

~/soft/MAVR/scripts/collaps_synonym_strings.py -i felis_catus.pep.sorted.description.no_isoform_versions.uniq \
                                               -k 1 -a 0 -o felis_catus.pep.collapsed_isoforms.description



