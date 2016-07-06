#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse

from Routines import SequenceRoutines, FileRoutines

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

args = parser.parse_args()

SequenceRoutines.make_region_bed_file_from_file(args.input, args.output, white_id_file=args.white_list_ids,
                                                black_id_file=args.black_list_ids,
                                                output_format="0-based", input_format=args.format)
