#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import EggNOGRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-d", "--dir_with_alignments", action="store", dest="dir_with_alignments", required=True,
                    help="Directory with alignments")
parser.add_argument("-o", "--output_dir", action="store", dest="output_dir", required=True,
                    help="Directory to write extracted proteins")

args = parser.parse_args()

EggNOGRoutines.extract_proteins_from_alignments(args.dir_with_alignments, args.output_dir)
