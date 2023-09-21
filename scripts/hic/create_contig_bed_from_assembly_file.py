#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
import pandas as pd


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input .assembly file")

parser.add_argument("-o", "--output", action="store", dest="output", default=sys.stdout,
                    help="Output five column (scaffold, start, end, strand, contig_id) .bed file with coordinates "
                         "of contigs in scaffolds. Default: stdout")

args = parser.parse_args()

with open(args.input, "r") as in_fd, open(args.output, "w") as out_fd:
    contig_df = []
    scaffold_list = []
    for line in in_fd:
        if line[0] == ">":
            contig_df.append(line[1:].strip().split())
        else:
            scaffold_list.append(line.strip().split())
    contig_df = pd.DataFrame.from_records(contig_df, columns=["contig_id", "cid", "length"], index=['cid'],)
    contig_df['length'] = contig_df['length'].astype(dtype='Int64')
    scaffold_index = 1

    out_fd.write("#scaffold_id\tstart\tend\tstrand\tcontig_id\n")
    for scaffold_components in scaffold_list:
        scaffold_id = "HiC_scaffold_{0}".format(scaffold_index)
        scaffold_shift = 0
        for c_____id in scaffold_components:
            if c_____id[0] == "-":
                strand = "-"
                cid = c_____id[1:]
            else:
                strand = "+"
                cid = c_____id
            contig_id = contig_df.loc[cid, "contig_id"]
            contig_len = contig_df.loc[cid, "length"]
            out_fd.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(scaffold_id,
                                                            scaffold_shift,
                                                            scaffold_shift + contig_len,
                                                            strand,
                                                            contig_id))
            scaffold_shift += contig_len

        scaffold_index += 1
