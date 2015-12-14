#!/usr/bin/env bash



~/soft/MAVR/scripts/annotation/extract_top_hits_from_target_gff_exonerate.py -i exonerate_cat_pep_splited_output.target.gff \
                                                                             -o exonerate_cat_pep_merged.rmdup \
                                                                             -d ~/data/genomes/cat/v8.0/pep/felis_catus.pep.sorted.description.uniq.ids

~/soft/augustus-3.2.1/scripts/exonerate2hints.pl --in=exonerate_cat_pep_merged.rmdup.target.top_hits.gff \
                                                 --out=exonerate_cat_pep_merged.rmdup.target.top_hits.hints \
                                                 --priority=20

~/soft/augustus-3.2.1/scripts/exonerate2hints.pl --in=exonerate_cat_pep_merged.rmdup.target.secondary_hits.gff \
                                                 --out=exonerate_cat_pep_merged.rmdup.target.secondary_hits.hints \
                                                 --priority=10

cat exonerate_cat_pep_merged.rmdup.target.top_hits.hints \
    exonerate_cat_pep_merged.rmdup.target.secondary_hits.hints > exonerate_cat_pep_merged.rmdup.target.merged.hints
