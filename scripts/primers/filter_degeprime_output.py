#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import PrimerRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with DEGEPRIME output")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file. Default: stdout")
parser.add_argument("-n", "--number_of_sequences", action="store", dest="num_of_seqs", required=True, type=int,
                    help="Total number of sequences")
parser.add_argument("-f", "--fraction", action="store", dest="fraction", type=float, default=0.9,
                    help="Minimum fraction of covered sequences to retain a primer. Default: 0.9")
parser.add_argument("-d", "--no_degeneration_dist", action="store", dest="no_deg_dist", type=int,
                    default=4,
                    help="Minimum number of 3' nucleotides without degeneration. Default: 4.")
parser.add_argument("-g", "--max_gc_content", action="store", dest="max_gc_content", type=float, default=0.75,
                    help="Maximum GC content of the primer. Default: 0.75")
parser.add_argument("-c", "--min_gc_content", action="store", dest="min_gc_content", type=float, default=0.4,
                    help="Minimum GC content of the primer. Default: 0.4")
parser.add_argument("-t", "--max_melting_temp_difference", action="store", dest="max_anneal_temp_diff",
                    type=float, default=6.0,
                    help="Maximum difference in melting temperature for degenerated primers. Default: 6.0 C")
parser.add_argument("-a", "--max_melting_temp", action="store", dest="max_melting_temp",
                    type=float, default=65.0,
                    help="Maximum melting temperature. Default: 65.0 C")
parser.add_argument("-e", "--min_melting_temp", action="store", dest="min_melting_temp",
                    type=float, default=50.0,
                    help="Minimum melting temperature. Default: 50.0 C")


args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

nucleotides = ["A", "C", "G", "T"]

low_coverage_counter = 0
degenerated_at_three_prime = 0
high_GC_count = 0
low_GC_count = 0
high_melting_temp_difference_count = 0
high_melting_temp_count = 0
low_melting_temp_count = 0

with open(args.input, "r") as in_fd:
    header = in_fd.readline()
    out_header = header.strip() + "\tMinGC\tMaxGC\tMim_Tm\tMaxTm\tTmDiff\n"
    out_fd.write(out_header)
    for line in in_fd:
        striped = line.strip()
        position, total_seq, unique_mers, entropy, degeneracy, matching_seq, primer_seq = striped.split("\t")
        if float(matching_seq)/float(args.num_of_seqs) < args.fraction:
            low_coverage_counter += 1
            continue
        for i in range(1, args.no_deg_dist + 1):
            #print i
            if primer_seq[-i] not in nucleotides:
                degenerated_at_three_prime += 1
                break
        else:
            primer_length = float(len(primer_seq))
            GC_AT_composition = PrimerRoutines.count_GC_AT_composition(primer_seq)

            max_GC_content = float(GC_AT_composition["GC"][1])/primer_length
            min_GC_content = float(GC_AT_composition["GC"][0])/primer_length

            if max_GC_content > args.max_gc_content:
                high_GC_count += 1
                continue
            if min_GC_content < args.min_gc_content:
                low_GC_count += 1
                continue

            melting_temperatures = PrimerRoutines.melting_temperature_formula2(primer_seq)

            melting_difference = melting_temperatures[1] - melting_temperatures[0]

            if melting_difference > args.max_anneal_temp_diff:
                high_melting_temp_difference_count += 1
                continue

            if melting_temperatures[1] > args.max_melting_temp:
                high_melting_temp_count += 1
                continue

            if melting_temperatures[0] < args.min_melting_temp:
                low_melting_temp_count += 1
                continue

            out_fd.write("%s\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n" % (striped, min_GC_content, max_GC_content,
                                                                 melting_temperatures[0], melting_temperatures[1],
                                                                 melting_difference))

out_fd.close()


