#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff",
                    help="Gff file with annotations to extract")
parser.add_argument("-o", "--output_gff", action="store", dest="output_gff",
                    help="Output gff file with extracted transcripts")
parser.add_argument("-d", "--ids_file", action="store", dest="ids_file",
                    help="File with ids of transcripts to extract")

args = parser.parse_args()

AnnotationsRoutines.extract_transcripts_by_ids(args.input_gff, args.ids_file, args.output_gff)
