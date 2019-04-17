#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.RepeatMasking import RepeatMasker

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with RepeatMasker output")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix for output files")
parser.add_argument("-e", "--expanded_set", action="store_true", dest="expanded_set", default=False,
                    help="Use expanded set of repeats. Includes 'Unknown', 'Other' and 'ARTEFACT' classes."
                         "Default: False")
parser.add_argument("-r", "--retain_unknown", action="store_true", dest="unknown", default=False,
                    help="Retain repeats from 'Unknown' class."
                         "Default: False")

args = parser.parse_args()

output_gff = "%s.gff" % args.output_prefix
repeat_classes = "%s.repeat_classes" % args.output_prefix
repeat_families = "%s.repeat_families" % args.output_prefix
extracted_gff = "%s.selected_repeat_classes.gff" % args.output_prefix

RepeatMasker.convert_rm_out_to_gff(args.input, output_gff, repeat_classes, repeat_families)
RepeatMasker.extract_repeats_used_for_gene_annotation(output_gff, extracted_gff, use_expanded_set=args.expanded_set,
                                                      retain_unknown_repeats=args.unknown)
