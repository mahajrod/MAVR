#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines, FileRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    type=lambda x: FileRoutines.make_list_of_path_to_files(x.split(",")),
                    help="Comma-separated list of genbank files/directories")
parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file")

args = parser.parse_args()

SequenceRoutines.extract_exon_lengths_from_genbank_file(args.input, args.output)



