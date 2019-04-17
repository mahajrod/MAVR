#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import AUGUSTUS


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input gff from AUGUSTUS")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-s", "--species_prefix", action="store", dest="species_prefix", required=True,
                    help="Species prefix for ids")
parser.add_argument("-d", "--number_of_digits_in_id", action="store", dest="number_of_digits_in_id", type=int,
                    default=8, help="Number of digits in ids. Default - 8")

args = parser.parse_args()

AUGUSTUS.replace_augustus_ids(args.input_gff, args.output_prefix, species_prefix=args.species_prefix,
                              number_of_digits_in_id=args.number_of_digits_in_id)
