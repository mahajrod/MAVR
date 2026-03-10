#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
from typing import OrderedDict

import pandas as pd
import argparse
from RouToolPa.Routines import AlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input filtered (no secondary alignments) alignment file with cigar string")
parser.add_argument("-f", "--format", action="store", dest="format", default="paf",
                    help="Format of the alignment file. Allowed: paf(default)")
parser.add_argument("-q", "--min_mapq", action="store", dest="min_mapq", default=20, type=int,
                    help="Minimum alignment quality. Default: 20")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

if args.format == "paf":
    aln_df = pd.read_csv(args.input, sep="\t", usecols=list(range(0,12)), header=None,
                         names=("qName", "qLen", "qStart", "qEnd", "qStrand",
                                "tName", "tLen", "tStart", "tEnd",
                                "matches", "alnLen", "MAPQ"))
    
    mapq_index = 11
    cigar_code_index_dict = {
                             "M": 12,
                             "I": 13,
                             "D": 14,
                             "N": 15,
                             "S": 16,
                             "H": 17,
                             "P": 18,
                             "=": 19,
                             "X": 20
                            }

    cigar_code_tuple = ("M", "I", "D", "N", "S", "H", "P", "=", "X")
    for cigar_code in cigar_code_tuple:
            aln_df[cigar_code] = 0

    line_index = 0
    with open(args.input, "r") as in_fd:
        for line in in_fd:
            line_list = line.strip().split()
            for element in line_list:
                if element[0:2] == "cg":
                    cigar_string = element.split(":")[-1]
                    cigar_code_dict = AlignmentRoutines.count_cigar_codes(AlignmentRoutines.parse_cigar(element.split(":")[-1]))
                    for cigar_code in cigar_code_dict:
                        aln_df.iloc[line_index, cigar_code_index_dict[cigar_code]] = cigar_code_dict[cigar_code]
                    break
            else:
                sys.stderr.write(f"WARNING!!! CIGAR string is absent in row (0-based) {line_index}\n !")

            line_index += 1

else:
    sys.stderr.write(f"ERROR!!! Unknown alignment format ({args.format})! Aborting...")
    exit(1)

aln_df.to_csv(f"{args.output_prefix}.cigar_stats.tab", sep="\t", header=True, index=False)

# filtering by mapping quality
raw_aln_count = len(aln_df)
aln_df = aln_df[aln_df[aln_df.columns[mapq_index]] >= args.min_mapq]
filtered_aln_count = len(aln_df)
aln_df.to_csv(f"{args.output_prefix}.cigar_stats.mapq{args.min_mapq}.tab", sep="\t", header=True, index=False)

count_dict = OrderedDict([(cigar_code, sum(aln_df[cigar_code])) for cigar_code in cigar_code_tuple])

cigar_match_code = None
# detect match code. It maybe either 'M' or '='
if (count_dict["M"] > 0) and (count_dict["="] > 0):
    sys.stderr.write(f"WARNING!!! CIGAR MATCH code is unclear. "
                     f"Both 'M' ({count_dict['M']}) and '=' ({count_dict['=']}) were detected in the file.")
elif (count_dict["M"] == 0) and (count_dict["="] == 0):
    sys.stderr.write(f"WARNING!!! CIGAR MATCH codes are absent in the file. "
                     f"No 'M' or '=' codes were detected in the file. ")
elif (count_dict["M"] > 0):
    cigar_match_code = "M"
elif (count_dict["="] > 0):
    cigar_match_code = "="

with open(f"{args.output_prefix}.all_stats.tab", "w") as stat_fd:

    stat_fd.write(f"Raw alignments\t{raw_aln_count}\n")
    stat_fd.write(f"MAPQ threshold\t{args.min_mapq}\n")
    stat_fd.write(f"Filtered alignments\t{filtered_aln_count}\n")
    stat_fd.write("CIGAR code counts:\n")

    for cigar_code in count_dict:
        stat_fd.write(f"{cigar_code}\t{count_dict[cigar_code]}\n")
    stat_fd.write("\n")

    if cigar_match_code is None:
        sys.stderr.write(f"WARNING!!! Skipping similarity calculations due to issues with the CIGAR match code")
    else:
        # base_similarity = matches / (matches + mismatches)
        base_similarity = float(count_dict[cigar_match_code]) / float(count_dict[cigar_match_code] + count_dict["X"])
        # similarity = matches / (matches + mismatches + insertions + deletions)
        similarity = float(count_dict[cigar_match_code]) / float(count_dict[cigar_match_code] + count_dict["X"] + count_dict["I"] + count_dict["D"])
        stat_fd.write(f"Base similarity(no indels): {base_similarity: .4f}\n")
        stat_fd.write(f"Similarity(including indels): {similarity: .4f}\n")