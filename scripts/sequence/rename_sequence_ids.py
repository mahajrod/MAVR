#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines.File import make_list_of_path_to_files
from RouToolPa.Routines import SequenceRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file_list", action="store", dest="input", required=True,
                    type=lambda s: make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of input files/directories with sequences")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with renamed sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-s", "--syn_file", action="store", dest="syn_file",
                    help="File with synonyms of ids")
parser.add_argument("-k", "--key_index", action="store", dest="key_index", type=int, default=0,
                    help="Key column in file with synonyms(0-based)")
parser.add_argument("-v", "--value_index", action="store", dest="value_index", type=int, default=1,
                    help="Value column in file with synonyms(0-based)")
parser.add_argument("-c", "--comments_prefix", action="store", dest="comments_prefix",
                    help="Prefix of comments in synonyms file")
parser.add_argument("-m", "--columns_separator", action="store", dest="separator", default="\t",
                    help="Column separator in file with synonyms")
parser.add_argument("-e", "--header", action="store_true", dest="header", default=False,
                    help="Header is present in synonyms file. Default - False")
parser.add_argument("-l", "--clear_description", action="store_true", dest="clear_description", default=False,
                    help="Clear description. Default - False")
parser.add_argument("-r", "--store_old_name_in_description", action="store_true", dest="store_old_name_in_description",
                    default=False,
                    help="Store old name in description. Default: False")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="parse",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db', 'index', 'parse'(default)")


args = parser.parse_args()

SequenceRoutines.rename_records_from_files(args.input, args.output, args.syn_file, format=args.format,
                                           header=args.header, separator=args.separator, key_index=args.key_index,
                                           value_index=args.value_index, clear_description=args.clear_description,
                                           store_old_name_in_description=args.store_old_name_in_description,
                                           comments_prefix=args.comments_prefix,  parsing_mode=args.parsing_mode)



