#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
import numpy as np
from Bio import AlignIO


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with alignment")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with substitutions. Default: stdout")
parser.add_argument("-r", "--reference_sequence_id", action="store", dest="ref_seq_id", required=True,
                    help="Id of reference sequence")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of alignment")
parser.add_argument("-g", "--gap_symbol", action="store", dest="gap_symbol", default="-",
                    help="Gap symbol used in alignment. Default: '-' ")
args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")
alignment = AlignIO.read(args.input, args.format)
number_of_seqs = len(alignment)

index_id_list = []
gaps_count_list = []
substitutions_list = []
for index in range(0, number_of_seqs):
    index_id_list.append(alignment[index].id)
    if alignment[index].id == args.ref_seq_id:
        reference_seq_index = index
    gaps_count_list.append(0)
    substitutions_list.append([])

alignment_array = np.array([list(record) for record in alignment], np.character, order="F")
length_of_alignment = np.size(alignment_array, 1)


for position in range(0, length_of_alignment):
    for seq_index in range(0, number_of_seqs):
        if alignment_array[seq_index, position] == args.gap_symbol:
            gaps_count_list[seq_index] += 1
        if seq_index == reference_seq_index:
            continue
        if alignment_array[seq_index, position] != alignment_array[reference_seq_index, position]:
            pos = position + 1 - gaps_count_list[reference_seq_index]

            if (alignment_array[reference_seq_index, position] != args.gap_symbol) and (alignment_array[seq_index, position] != args.gap_symbol):
                substitution = "%s%i%s" % (alignment_array[reference_seq_index, position], pos,
                                           alignment_array[seq_index, position])
            elif alignment_array[reference_seq_index, position] == args.gap_symbol:
                substitution = "Ins%i%s" % (pos, alignment_array[seq_index, position])
            else:
                substitution = "Del%s%i" % (alignment_array[reference_seq_index, position], pos)
            substitutions_list[seq_index].append(substitution)

for index in range(0, number_of_seqs):
    if index == reference_seq_index:
        continue
    out_fd.write("%s\t%s\n" % (index_id_list[index],
                               ",".join(substitutions_list[index]) if substitutions_list[index] else "."))
if args.output != "stdout":
    out_fd.close()