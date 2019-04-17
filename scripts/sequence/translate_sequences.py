#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import FileRoutines, SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=lambda s: FileRoutines.make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of input files/directories with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output", required=True,
                    help="Output file with translated sequences")
parser.add_argument("-g", "--genetic_code", action="store", dest="code", type=int, default=1,
                    help="Genetic code to use. Set by number of ncbi code. Default - 1 (universal)")
parser.add_argument("-s", "--translate_to_stop", action="store_true", dest="translate_to_stop",
                    help="Translate to first in-frame stop codon. Default - False(translate whole sequence)")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="parse",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db', 'index', 'parse'(default)")

args = parser.parse_args()


SequenceRoutines.translate_sequences_from_file(args.input, args.output, format="fasta", id_expression=None,
                                               genetic_code_table=args.code, translate_to_stop=True,
                                               mode=args.parsing_mode, index_file="tmp.idx")
