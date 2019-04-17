#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Collections.General import IdList, SynDict

# from RouToolPa.Routines.File import read_synonyms_dict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with families")
parser.add_argument("-d", "--id_file", action="store", dest="id_file", default=None,
                    help="File with ids of families. If absent genes from all families will be extracted(default).")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file")
parser.add_argument("-s", "--separate_families", action="store_true", dest="separate_families",
                    help="Separate families to different files. If set option -o/--output_file is ignored")
args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
families = SynDict()
families.read(args.input, separator="\t", split_values=True, values_separator=",")
if args.id_file:
    id_list = IdList()
    id_list = id_list.read(args.id_file)


if args.separate_families:
    for fam_id in id_list if args.id_file else families:
        with open("%s.ids" % fam_id, "w") as fam_fd:
            for gene_id in families[fam_id]:
                fam_fd.write(gene_id + "\n")
else:
    with open(args.output, "w") as out_fd:
        for fam_id in id_list if args.id_file else families:
            for gene_id in families[fam_id]:
                out_fd.write(gene_id + "\n")
if args.output != "stdout":
    out_fd.close()
