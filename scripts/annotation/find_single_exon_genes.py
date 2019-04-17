#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from Bio import SeqIO
from BCBio import GFF
from RouToolPa.Collections.General import IdList



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input", required=True,
                    help="Input .gff file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output", default="stdout",
                    help="Output file with single exon genes. Default: stdout")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
annotations_dict = SeqIO.to_dict(GFF.parse(open(args.input)))
single_gene_id_list = IdList()

for record in annotations_dict:
    for feature in annotations_dict[record].features:
        #print feature.id
        if feature.type != "gene":
            continue
        for subfeature in feature.sub_features:
            if subfeature.type != "mRNA":
                continue
            exon_number = 0
            for mRNA_subfeature in subfeature.sub_features:
                if mRNA_subfeature.type == "exon":
                    exon_number += 1
            if exon_number == 1:
                single_gene_id_list.append(feature.id)

single_gene_id_list.write(out_fd, close_after_if_file_object=True)
"""
sequence_groups_id = SynDict()
sequence_groups_id.read(args.id_file, split_values=True)
#print("Parsing %s..." % args.input_file)
sequence_dict = SeqIO.index_db(tmp_index_file, args.input, format=args.format)
for group in sequence_groups_id:
    SeqIO.write(record_by_id_generator(sequence_dict, sequence_groups_id[group]),
                "%s%s.%s" % (args.output, group, args.extension), format=args.format)
"""





