#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with sequences")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output file")
parser.add_argument("-r", "--regions", action="store", dest="regions", required=True,
                    help="File with regions. Default column indexes are adjusted for BED files")
parser.add_argument("-m", "--min_length", action="store", dest="min_length",
                    type=int, default=1,
                    help="Minimum length of splited region to retain")

parser.add_argument("-c", "--scaffold_column_id", action="store", dest="scaffold_column_id",
                    type=int, default=0,
                    help="0-based index of column with scaffold id. Default: 0")
parser.add_argument("-s", "--start_column_id", action="store", dest="start_column_id",
                    type=int, default=1,
                    help="0-based index of column with feature start. Default: 1")
parser.add_argument("-e", "--end_column_id", action="store", dest="end_column_id",
                    type=int, default=2,
                    help="0-based index of column with feature end. Default: 2")
parser.add_argument("-n", "--coordinates_type", action="store", dest="coordinates_type",
                    default="1-based",
                    help="Type of coordinates. Allowed: 0-based, 1-based(default)")

args = parser.parse_args()

SequenceRoutines.split_sequence_by_regions_from_file(args.input,
                                                     args.regions,
                                                     args.output_prefix,
                                                     retain_description=False,
                                                     min_length=args.min_length,
                                                     parsing_mode="parse",
                                                     scaffold_column_index=args.scaffold_column_id,
                                                     start_column_index=args.start_column_id,
                                                     end_column_index=args.end_column_id,
                                                     coordinates_type=args.coordinates_type,
                                                     input_separator="\t",
                                                     sequence_format="fasta")
