#!/usr/bin/env bash

#get longest proteins per genome for eshark
grep -P "^>" eshark.pep.fa | sed 's/>//;s/len=//' | awk '{printf "%s\t%s\t%s\n", $4, $1, $2}' | sort -k1,1 -k3,3nr > eshark.pep.len.tab

awk '{if ($1 != PREV) {print $0; PREV=$1}}' eshark.pep.len.tab  > eshark.pep.len.longest_pep.tab

awk '{print $2}' eshark.pep.len.longest_pep.tab > eshark.pep.len.longest_pep.ids




#get longest proteins for whale shark

grep -P '\tCDS\t' /home/mahajrod/Genetics/Projects/white_shark/species/rhincodon_typus/data/annotation/GCF_001642345.1_ASM164234v2_genomic.gff | grep 'protein_id=' > GCF_001642345.1_ASM164234v2_genomic.cds.gff

sed 's/.*gene=\([A-Za-z0-9\.\-\_]\+\).*protein_id=\([A-Za-z0-9\.\-\_]\+\).*/\1\t\2/' GCF_001642345.1_ASM164234v2_genomic.cds.gff | sort | uniq > GCF_001642345.1_ASM164234v2_genomic.gene_to_protein.accordance