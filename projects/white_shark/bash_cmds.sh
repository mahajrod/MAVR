#!/usr/bin/env bash

#get longest proteins per genome for eshark
grep -P "^>" eshark.pep.fa | sed 's/>//;s/len=//' | awk '{printf "%s\t%s\t%s\n", $4, $1, $2}' | sort -k1,1 -k3,3nr > eshark.pep.len.tab

awk '{if ($1 != PREV) {print $0; PREV=$1}}' eshark.pep.len.tab  > eshark.pep.len.longest_pep.tab

awk '{print $2}' eshark.pep.len.longest_pep.tab > eshark.pep.len.longest_pep.ids

#get longest proteins for whale shark

grep -P '\tCDS\t' /home/mahajrod/Genetics/Projects/white_shark/species/rhincodon_typus/data/annotation/GCF_001642345.1_ASM164234v2_genomic.gff | grep 'protein_id=' > GCF_001642345.1_ASM164234v2_genomic.cds.gff

sed 's/.*gene=\([A-Za-z0-9\.\-\_]\+\).*protein_id=\([A-Za-z0-9\.\-\_]\+\).*/\1\t\2/' GCF_001642345.1_ASM164234v2_genomic.cds.gff | sort | uniq > GCF_001642345.1_ASM164234v2_genomic.gene_to_protein.accordance

#Astyanax mexicanus
~/Soft/MAVR/scripts/annotation/ensembl/get_gene_transcript_protein_from_ensembl_pep_fasta.py -i Astyanax_mexicanus.AstMex102.pep.all.fa -o astyanax_mexicanus.gene.trascript.protein.tsv

~/Soft/MAVR/scripts/file/extract_by_column_value.py -i astyanax_mexicanus.gene.trascript.protein.tsv -c 2 -f Astyanax_mexicanus.longest_pep.ids -o astyanax_mexicanus.gene.trascript.protein.longest_pep

awk -F'\t' '{printf "%s\t%s\n", $2, $3}' astyanax_mexicanus.gene.trascript.protein.longest_pep.tsv > astyanax_mexicanus.cds_to_pep.accordance

awk -F'\t' '{print $2}' astyanax_mexicanus.gene.trascript.protein.longest_pep.tsv > astyanax_mexicanus.longest_cds.ids

~/Soft/MAVR/scripts/annotation/trim_cds_and_remove_terminal_stop_codons.py -i Astyanax_mexicanus.AstMex102.cds.all.fa  -o Astyanax_mexicanus.AstMex102.cds.all.trimmed.fa
~/Soft/MAVR/scripts/sequence/extract_sequences_by_ids.py -i Astyanax_mexicanus.AstMex102.cds.all.trimmed.fa -o Astyanax_mexicanus.longest_cds.cds -d astyanax_mexicanus.longest_cds.ids

#Poecilia_formosa

~/Soft/MAVR/scripts/annotation/ensembl/get_gene_transcript_protein_from_ensembl_pep_fasta.py -i Poecilia_formosa.PoeFor_5.1.2.pep.all.fa -o poecilia_formosa.gene.trascript.protein.tsv

~/Soft/MAVR/scripts/file/extract_by_column_value.py -i poecilia_formosa.gene.trascript.protein.tsv -c 2 -f Poecilia_formosa.longest_pep.ids -o poecilia_formosa.gene.trascript.protein.longest_pep

awk -F'\t' '{printf "%s\t%s\n", $2, $3}' poecilia_formosa.gene.trascript.protein.longest_pep.tsv > poecilia_formosa.cds_to_pep.accordance

awk -F'\t' '{print $2}' poecilia_formosa.gene.trascript.protein.longest_pep.tsv > poecilia_formosa.longest_cds.ids

~/Soft/MAVR/scripts/annotation/trim_cds_and_remove_terminal_stop_codons.py -i Poecilia_formosa.PoeFor_5.1.2.cds.all.fa  -o Poecilia_formosa.PoeFor_5.1.2.cds.all.trimmed.fa
~/Soft/MAVR/scripts/sequence/extract_sequences_by_ids.py -i Poecilia_formosa.PoeFor_5.1.2.cds.all.trimmed.fa -o Poecilia_formosa.longest_cds.cds -d poecilia_formosa.longest_cds.ids

#Lepisosteus_oculatus

~/Soft/MAVR/scripts/annotation/ensembl/get_gene_transcript_protein_from_ensembl_pep_fasta.py -i Lepisosteus_oculatus.LepOcu1.pep.all.fa -o lepisosteus_oculatus.gene.trascript.protein.tsv

~/Soft/MAVR/scripts/file/extract_by_column_value.py -i lepisosteus_oculatus.gene.trascript.protein.tsv -c 2 -f Lepisosteus_oculatus.longest_pep.ids -o lepisosteus_oculatus.gene.trascript.protein.longest_pep

awk -F'\t' '{printf "%s\t%s\n", $2, $3}' lepisosteus_oculatus.gene.trascript.protein.longest_pep.tsv > lepisosteus_oculatus.cds_to_pep.accordance

awk -F'\t' '{print $2}' lepisosteus_oculatus.gene.trascript.protein.longest_pep.tsv > lepisosteus_oculatus.longest_cds.ids

~/Soft/MAVR/scripts/annotation/trim_cds_and_remove_terminal_stop_codons.py -i Lepisosteus_oculatus.LepOcu1.cds.all.fa  -o Lepisosteus_oculatus.LepOcu1.cds.all.trimmed.fa
~/Soft/MAVR/scripts/sequence/extract_sequences_by_ids.py -i Lepisosteus_oculatus.LepOcu1.cds.all.trimmed.fa -o Lepisosteus_oculatus.longest_cds.cds -d lepisosteus_oculatus.longest_cds.ids

#callorhinchus_milii
grep -P "^>" eshark.pep.fa | sed 's/>//;s/\([A-Za-z0-9\._-]\+\)\s.*\s\([A-Za-z0-9\._-]\+\)\s\([A-Za-z0-9\._-]\+\)/\3\t\2\t\1/' > callorhinchus_milii.gene.trascript.protein.tsv
awk -F'\t' '{printf "%s\t%s\n", $2, $3}' callorhinchus_milii.gene.trascript.protein.tsv > callorhinchus_milii.transcript_to_pep.accordance
~/Soft/MAVR/scripts/sequence/find_cds_coordinates_in_transcript_by_pep.py -t eshark_transcripts.fa  -p eshark.pep.fa  -o callorhinchus_milii.cds.coords -c callorhinchus_milii.transcript_to_pep.accordance
~/Soft/MAVR/scripts/transcriptome/annotation/get_cds_by_bed_from_transcript_file.py  -i eshark_transcripts.fa -o callorhinchus_milii.cds  -c callorhinchus_milii.cds.coords

~/Soft/MAVR/scripts/annotation/trim_cds_and_remove_terminal_stop_codons.py -i callorhinchus_milii.cds  -o callorhinchus_milii.trimmed.cds

~/Soft/MAVR/scripts/file/extract_by_column_value.py -i callorhinchus_milii.gene.trascript.protein.tsv -c 2 -f callorhinchus_milii.longest_pep.ids -o callorhinchus_milii.gene.trascript.protein.longest_pep

awk -F'\t' '{print $2}' callorhinchus_milii.gene.trascript.protein.longest_pep.tsv > callorhinchus_milii.longest_cds.ids

~/Soft/MAVR/scripts/sequence/extract_sequences_by_ids.py -i callorhinchus_milii.trimmed.cds -o Callorhinchus_milii.longest_cds.cds -d callorhinchus_milii.longest_cds.ids

#rhincodon_typus

~/Soft/MAVR/scripts/sequence/get_sequence_ids.py -i rhincodon_typus.longest.pep -o rhincodon_typus.longest.pep.ids
~/Soft/MAVR/scripts/file/extract_by_column_value.py -i rhincodon_typus.cds.pep_to_transcript.accordance -c 0 -f rhincodon_typus.longest.pep.ids -o rhincodon_typus.pep_to_transcript.longest_pep.accordance

awk -F'\t' '{printf "%s\t%s\n", $2, $1}' rhincodon_typus.pep_to_transcript.longest_pep.tsv.accordance.accordance > rhincodon_typus.transcript_to_pep.longest_pep.accordance
awk -F'\t' '{print $2}' rhincodon_typus.transcript_to_pep.longest_pep.accordance > rhincodon_typus.longest_cds.ids
~/Soft/MAVR/scripts/sequence/extract_sequences_by_ids.py -i rhincodon_typus.cds.cds -o Rhincodon_typus.longest_cds.cds -d rhincodon_typus.longest_cds.ids



