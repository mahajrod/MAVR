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
                    help="Output file with sequences")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input and output files. Allowed formats genbank, fasta(default)")
parser.add_argument("-d", "--id_file", action="store", dest="id_file",
                    help="File with ids or comma-separated list of ids of sequences to extract")
parser.add_argument("-c", "--id_column", action="store", dest="id_column", type=int, default=0,
                    help="Number of column with ids in id file (0-based). Default: 0")
parser.add_argument("-e", "--extraction_mode", action="store", dest="coincidence_mode", default="exact",
                    help="Coincidence mode for id: exact(full, default), partial")
parser.add_argument("-a", "--allow_multiple_coincidence_report", action="store_true",
                    dest="allow_multiple_coincidence_report", default=False,
                    help="Allow multiple coincidence report of sequences for partial coincidence mode."
                         "By default an error is raised")
parser.add_argument("-s", "--syn_file", action="store", dest="syn_file",
                    help="File with synonyms of ids to use. Default - not set")
parser.add_argument("-r", "--invert_match", action="store_true", dest="invert_match",
                    help="Invert match, i. e. remove sequences. Default - not set")

args = parser.parse_args()

SequenceRoutines.extract_sequence_by_ids(args.input, args.id_file, args.output, format=args.format, verbose=True,
                                         id_column_number=args.id_column, coincidence_mode=args.coincidence_mode,
                                         allow_multiple_coincidence_report=args.allow_multiple_coincidence_report,
                                         syn_file=args.syn_file, parsing_mode="parse", index_file="tmp.idx",
                                         invert_match=args.invert_match)
