#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Routines import AnnotationsRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-g", "--gff", action="store", dest="gff", required=True,
                    help="Gff file")

parser.add_argument("-o", "--output", action="store", dest="output",
                    help="Output file with ids. Default: stdout")

args = parser.parse_args()

if args.output is None:
    args.output = sys.stdout

AnnotationsRoutines.get_scaffold_ids_from_gff(args.gff, out_file=args.output)
