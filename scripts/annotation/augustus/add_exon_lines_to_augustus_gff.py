#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import AUGUSTUS

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input AUGUSTUS GFF file")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff", required=True,
                    help="Output GFF with exon entries")
parser.add_argument("-e", "--exon_id_prefix", action="store", dest="exon_id_prefix", default="EXON",
                    help="Prefix of exon id. Default: EXON")
parser.add_argument("-n", "--id_digit_num", action="store", dest="id_digit_num", default=8, type=int,
                    help="Number of digits in exon id. Default: 8")


args = parser.parse_args()

AUGUSTUS.add_exon_lines_to_augustus_gff(args.input_gff, args.output_gff, number_of_digits_in_id=args.id_digit_num,
                                        exon_id_prefix=args.exon_id_prefix,
                                        new_exon_numering=False)
