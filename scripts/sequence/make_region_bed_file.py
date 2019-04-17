#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines, FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    type=lambda s: FileRoutines.make_list_of_path_to_files(s.split(",")),
                    help="Comma-separated list of files/directories with sequences")
parser.add_argument("-o", "--output_bed_file", action="store", dest="output", required=True,
                    help="Output bed file")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input files. Allowed formats genbank, fasta(default)")
parser.add_argument("-w", "--white_list_ids", action="store", dest="white_list_ids",
                    help="File with ids of regions from white list")
parser.add_argument("-b", "--black_list_ids", action="store", dest="black_list_ids",
                    help="File with ids of regions from black list")
parser.add_argument("-e", "--bed_format", action="store", dest="bed_format", default="0-based",
                    help="Format of output bed format. Allowed: 0-based(default), 1-based")
parser.add_argument("-m", "--min_len", action="store", dest="min_len",
                    help="Minimum length of sequence to count. Default: not set")
parser.add_argument("-x", "--max_len", action="store", dest="max_len",
                    help="Maximum length of sequence to count. Default: not set")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="index_db",
                    help="Parsing mode for input sequence file. "
                         "Possible variants: 'index_db'(default), 'index', 'parse'")

args = parser.parse_args()

SequenceRoutines.make_region_bed_file_from_file(args.input, args.output, white_id_file=args.white_list_ids,
                                                black_id_file=args.black_list_ids,
                                                output_format=args.bed_format, input_format=args.format,
                                                min_len=args.min_len, max_len=args.max_len,
                                                parsing_mode=args.parsing_mode,
                                                index_file="tmp.idx", retain_index=False)
