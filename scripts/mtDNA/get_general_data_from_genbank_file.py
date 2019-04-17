#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import os
from Bio import SeqIO
from RouToolPa.Routines import SequenceRoutines #  record_by_id_generator, find_gaps
from GeneSynonyms.MitoGenesSynonyms import mithochondrioal_genes_syndict



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input genbank file with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="out_file",
                    help="Output file")


args = parser.parse_args()


tmp_index_file = "temp.idx"

print("Parsing %s..." % args.input_file)
sequence_dict = SeqIO.index_db(tmp_index_file, args.input_file, format="genbank")
gaps_dict = SequenceRoutines.find_gaps(sequence_dict)

with open(args.out_file, "w") as out_fd:
    out_fd.write("#sequence_id\tspecies\tlength\tlineage\treferences\tgenes\trRNA\n")
    for record_id in sequence_dict:
        species = sequence_dict[record_id].annotations["organism"]
        lineage = sequence_dict[record_id].annotations["taxonomy"]
        length = len(sequence_dict[record_id].seq)
        references = ",".join(map(lambda x: x.title, sequence_dict[record_id].annotations['references']))

        protein_genes_list = []
        rRNA_genes_list = []
        for feature in sequence_dict[record_id].features:
            if feature.type == "gene":
                if "gene" in feature.qualifiers:
                    gene_name = feature.qualifiers["gene"][0]
                elif "product" in feature.qualifiers:
                    gene_name = feature.qualifiers["product"][0]
                else:
                    gene_name = "unnamed_gene"

                for mt_gene_name in mithochondrioal_genes_syndict:
                    if gene_name in mithochondrioal_genes_syndict[mt_gene_name]:
                        gene_name = mt_gene_name
                        break

                protein_genes_list.append(gene_name)
            elif feature.type == "rRNA":
                if "gene" in feature.qualifiers:
                    rRNA_name = feature.qualifiers["gene"][0]
                elif "product" in feature.qualifiers:
                    rRNA_name = feature.qualifiers["product"][0]
                else:
                    rRNA_name = "unnamed_rRNA"
                for mt_gene_name in mithochondrioal_genes_syndict:
                    if rRNA_name in mithochondrioal_genes_syndict[mt_gene_name]:
                        rRNA_name = mt_gene_name
                        break
                rRNA_genes_list.append(rRNA_name)

        out_fd.write("%s\t%s\t%i\t%s\t%s\t%s\t%s\n" % (record_id, species, length, ",".join(lineage), references,
                                                        ",".join(protein_genes_list), ",".join(rRNA_genes_list)))

awk_string = "awk -F'\\t' 'NR==1 {}; NR>1 {printf \"%%s\\t%%s\\n\", $2, $3}' %s | sort -t $'\\t' -k 1 -k2n > %s" % (args.out_file,
                                                                                         args.out_file + ".len")
os.remove(tmp_index_file)
