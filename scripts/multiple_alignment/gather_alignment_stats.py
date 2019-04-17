#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Parsers.MultipleAlignment import MultipleAlignmentStatCollection

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", type=lambda x: x.split(","), required=True,
                    help="Comma-separated list of files or directory with files containing sequences with alignment")
"""
parser.add_argument("-o", "--output_directory", action="store", dest="output", type=FileRoutines.check_path,
                    help="Output directory")
"""
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Verbose parsing. Default: False")


"""
parser.add_argument("-m", "--mode", action="store", dest="mode", default="globalpair",
                    help="Alignment mode. Default: 'globalpair'. Allowed: globalpair, localpair, genafpair")
parser.add_argument("-q", "--quiet", action="store_true", dest="quiet",
                    help="Quiet output")
parser.add_argument("-x", "--maxiterate", action="store", dest="maxiterate", type=float,
                    help="Maximum number of iterations")
parser.add_argument("-f", "--offset", action="store", dest="offset", type=float,
                    help="Offset (works like gap extension penalty)")
parser.add_argument("-g", "--gap_open_penalty", action="store", dest="gap_open_penalty", type=float,
                    help="Gap open penalty")
"""
args = parser.parse_args()

general_stat_file = "%s.general.stats" % args.output_prefix

stat_collection = MultipleAlignmentStatCollection(input_file=args.input, verbose=args.verbose)
stat_collection.write_general_stats(general_stat_file)
stat_collection.write_stats(args.output_prefix)
