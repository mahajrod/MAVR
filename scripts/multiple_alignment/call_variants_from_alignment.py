#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import MultipleAlignmentRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with alignment")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-r", "--reference_sequence_id", action="store", dest="reference_sequence_id",
                    required=True,
                    help="Reference sequence id")
parser.add_argument("-g", "--gap_symbol", action="store", dest="gap_symbol", default="-",
                    help="Gap symbol. Default - '-'")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Alignment file format. Default: fasta")
parser.add_argument("-a", "--align_variants", action="store_true", dest="align_variants",
                    help="Align variants between species by coordinate. Default: False")
parser.add_argument("-t", "--target_sequence_id", action="store", dest="target_sequence_id",
                    help="Target sequence id. Variants specific for this sequence will be "
                         "extracted into separated file. Default: not set")
args = parser.parse_args()

MultipleAlignmentRoutines.call_variants_from_multiple_alignment_from_file(args.input,
                                                                          args.output_prefix,
                                                                          args.reference_sequence_id,
                                                                          gap_symbol=args.gap_symbol,
                                                                          verbose=True,
                                                                          format="fasta",
                                                                          align_variants=args.align_variants,
                                                                          output_type="hgvs",
                                                                          variant_separator=",",
                                                                          target_sequence_id=args.target_sequence_id,
                                                                          absent_symbol="")
