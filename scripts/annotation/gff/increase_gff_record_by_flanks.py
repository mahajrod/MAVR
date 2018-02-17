#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse

from Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff", required=True,
                    help="Input .gff file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-f", "--fasta", action="store", dest="fasta", required=True,
                    help="Fasta file with sequences")
parser.add_argument("-d", "--id_description_entry", action="store",dest="id_description_entry", default="ID",
                    help="Key of feature_id entry ion description field. Default: ID ")
parser.add_argument("-c", "--coords_description_entry", action="store", dest="coords_description_entry",
                    default="core_seq_coords",
                    help="Key for description entry with coordinates of core sequence in new feature")

parser.add_argument("-l", "--left_flank_len", action="store", dest="left_flank_len", type=int, default=300,
                    help="Length of left flank. Default: 300")
parser.add_argument("-r", "--right_right_len", action="store", dest="right_flank_len", type=int, default=300,
                    help="Length of right flank. Default: 300")


args = parser.parse_args()

AnnotationsRoutines.increase_gff_record_by_flanks(args.input_gff, args.output_prefix,
                                                  args.left_flank_len, args.right_flank_len, args.fasta,
                                                  coords_description_entry=args.coords_description_entry,
                                                  id_description_entry=args.id_description_entry)
