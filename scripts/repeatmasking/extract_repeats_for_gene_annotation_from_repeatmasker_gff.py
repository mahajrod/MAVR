#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.RepeatMasking import RepeatMasker

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input .gff file with RepeatMasker output")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output gff with selected repeats")
parser.add_argument("-e", "--expanded_set", action="store_true", dest="expanded_set", default=False,
                    help="Use expanded set of repeats. Includes 'Unknown', 'Other' and 'ARTEFACT' classes."
                         "Default: False")
parser.add_argument("-r", "--retain_unknown", action="store_true", dest="unknown", default=False,
                    help="Retain repeats from 'Unknown' class."
                         "Default: False")

args = parser.parse_args()

RepeatMasker.extract_repeats_used_for_gene_annotation(args.input, args.output, use_expanded_set=args.expanded_set,
                                                      retain_unknown_repeats=args.unknown)
