#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
import os
from Bio import SeqIO
from RouToolPa.Routines.Sequence import find_gaps


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input genbank file with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="out_file",
                    help="Output file")


args = parser.parse_args()


tmp_index_file = "temp.idx"

print("Parsing %s..." % args.input_file)
sequence_dict = SeqIO.index_db(tmp_index_file, args.input_file, format="genbank")
gaps_dict = find_gaps(sequence_dict)

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

                protein_genes_list.append(gene_name)
            elif feature.type == "rRNA":
                if "gene" in feature.qualifiers:
                    rRNA_name = feature.qualifiers["gene"][0]
                elif "product" in feature.qualifiers:
                    rRNA_name = feature.qualifiers["product"][0]
                else:
                    rRNA_name = "unnamed_rRNA"

                rRNA_genes_list.append(rRNA_name)

        out_fd.write("%s\t%s\t%i\t%s\t%s\t%s\t%s\n" % (record_id, species, length, ",".join(lineage), references,
                                                        ",".join(protein_genes_list), ",".join(rRNA_genes_list)))

os.remove(tmp_index_file)



