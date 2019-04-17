#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Comma-separated list of files or directory with files containing sequences with alignment")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with dN, dS, W")
parser.add_argument("-r", "--ref_seq_id", action="store", dest="ref_seq_id", required=True,
                    help="Id of sequence to be used as reference")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Alignment format. Default: fasta")
parser.add_argument("-g", "--genetic_code_table", action="store", dest="genetic_code_table", type=int,
                    default=1,
                    help="Genetic code table number")
parser.add_argument("-a", "--gap_symbol_list", action="store", dest="gap_symbol_list", type=lambda s: s.split(","),
                    default=["-"],
                    help="Comma-separated list of gap symbols. Default: '-'")
parser.add_argument("-b", "--use_ambiguous_table", action="store_true", dest="use_ambiguous_table", default=False,
                    help="Use ambiguous codon table. Default:False")

args = parser.parse_args()


MultipleAlignmentRoutines.count_dNdS_by_reference_seq_in_codon_alignment_from_file(args.input,
                                                                                   args.ref_seq_id,
                                                                                   genetic_code_table=args.genetic_code_table,
                                                                                   gap_symbol_list=args.gap_symbol_list,
                                                                                   use_ambigious_table=args.ambigious_table,
                                                                                   output_file=args.output,
                                                                                   format="fasta")
