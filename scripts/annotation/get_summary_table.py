#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from collections import OrderedDict
from Bio.SeqUtils import seq1
from RouToolPa.Routines.File import make_list_of_path_to_files, split_filename
from RouToolPa.Collections.General import TwoLvlDict, SynDict



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", type=lambda s: s.split(","),
                    help="Comma-separated list of directories/files with SNPeff output extracted from vcf file")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file. Default: stdout")
parser.add_argument("-s", "--suffix_to_remove", action="store", dest="suffix_to_remove",
                    help="Suffix to remove from filename(not count extension)")
parser.add_argument("-g", "--gene_alias_file", action="store", dest="gene_alias_file",
                    help="File with aliases of genes (tab-separated, one alias per gene, no header)")
parser.add_argument("-w", "--write_dir_path", action="store_true", dest="write_dir_path",
                    help="write directory name(if directory is source of vcf files) in output file. Default: false")
parser.add_argument("-e", "--write_ext", action="store_true", dest="write_ext",
                    help="write extensions of vcf files in output file. Default: false")
parser.add_argument("-r", "--remove_nucleotide_substitutions", action="store_true", dest="rem_nuc_sub",
                    help="Remove nucleotide substitutions from output(preserve only AA substitutions)")
parser.add_argument("-c", "--convert_aa_to_single_letter", action="store_true", dest="convert_to_single_letter",
                    help="Convert aminoacids to single letters")

args = parser.parse_args()

args.input = make_list_of_path_to_files(args.input)

gene_alias_dict = SynDict()
if args.gene_alias_file:
    gene_alias_dict.read(args.gene_alias_file, split_values=False)
out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

summary_dict = TwoLvlDict()
for filename in args.input:
    directory, prefix, extension = split_filename(filename)

    if args.write_dir_path and args.write_ext:
        name = filename
    elif args.write_dir_path:
        name = (directory + prefix) if directory else prefix
    elif args.write_ext:
        name = prefix + extension
    else:
        name = prefix
        if args.suffix_to_remove in name:
            name = name.replace(args.suffix_to_remove, "")
    summary_dict[name] = OrderedDict()
    with open(filename, "r") as file_fd:
        file_fd.readline()
        for line in file_fd:
            tmp = line.strip().split("\t")
            gene_name = tmp[10]
            substitution = tmp[8]
            if substitution == ".":  # skip substitutions not in CDS
                continue
            if gene_alias_dict:
                if gene_name in gene_alias_dict:
                    gene_name = gene_alias_dict[gene_name]
            if args.rem_nuc_sub:
                substitution = substitution.split("/")[0][2:]
                if args.convert_to_single_letter:
                    ref_aa = seq1(substitution[:3])
                    try:
                        if substitution[-1] == "*":
                            alt_aa = "*"
                            pos = substitution[3:-1]
                        else:
                            alt_aa = seq1(substitution[-3:])
                            pos = substitution[3:-3]
                        substitution = "%s%s%s" % (ref_aa, pos, alt_aa)
                    except:
                        print(substitution, "aaa", filename, gene_name)
            if gene_name not in summary_dict[name]:
                summary_dict[name][gene_name] = [substitution]
            else:
                summary_dict[name][gene_name].append(substitution)

summary_dict.write(out_fd, absent_symbol=".")
if args.output != "stdout":
    out_fd.close()

