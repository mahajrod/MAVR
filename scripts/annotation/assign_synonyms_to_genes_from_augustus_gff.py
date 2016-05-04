#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from Tools.Annotation import AUGUSTUS

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input gff from AUGUSTUS")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file with synonyms")
parser.add_argument("-p", "--id_prefix", action="store", dest="id_prefix", required=True,
                    help="Prefix of id")
parser.add_argument("-n", "--number_of_digits_in_number", action="store", dest="number_of_digits_in_number", type=int,
                    default=6, help="Number of digits in id")

args = parser.parse_args()


AUGUSTUS.assign_synonyms_to_genes_from_augustus_gff(args.input_gff, args.output, args.id_prefix,
                                                    number_of_digits_in_number=args.number_of_digits_in_number)
