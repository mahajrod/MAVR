#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-t", "--transcript_file", action="store", dest="transcript_file", required=True,
                    help="Input file with transcript sequences")
parser.add_argument("-p", "--pep_file", action="store", dest="pep_file", required=True,
                    help="Input file with protein sequences")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output_files")
parser.add_argument("-c", "--correspondence_file", action="store", dest="correspondence_file",
                    help="File with transcript to protein correspondence")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input files. Allowed: fasta, genbank. Default: fasta")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Print warning if no protein was found for CDS")
parser.add_argument("-m", "--parsing_mode", action="store", dest="parsing_mode", default="parse",
                    help="Parsing mode of sequence files. Allowed: parse, index, index_db."
                         "Default: parse")
parser.add_argument("-g", "--genetic_code_table", action="store", dest="genetic_code_table", default=1, type=int,
                    help="Genetic code to use for translation of transcript. "
                         "Allowed: table number from http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"
                         "Default: 1(The standard code)")

args = parser.parse_args()

SequenceRoutines.find_cds_coordinates_in_transcript_by_pep_from_file(args.transcript_file, args.pep_file,
                                                                     args.correspondence_file, args.output_prefix,
                                                                     parsing_mode=args.parsing_mode,
                                                                     verbose=args.verbose,
                                                                     format=args.format, transcript_index_file=None,
                                                                     protein_index_file=None,
                                                                     genetic_code_table=args.genetic_code_table)
