#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import AnnotationsRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-c", "--correspondence_file", action="store", dest="correspondence_file", required=True,
                    help="File with correspondence of transcripts to genes")
parser.add_argument("-l", "--length_file", action="store", dest="length_file", required=True,
                    help="Length file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")


args = parser.parse_args()

AnnotationsRoutines.add_length_to_accordance_file(args.correspondence_file, args.length_file, args.output_prefix)
