#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import AnnotationsRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--gtf_file", action="store", dest="input", required=True,
                    help="Input gtf file")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output accordance file")

args = parser.parse_args()

AnnotationsRoutines.get_transcript_to_pep_accordance_from_gtf(args.input, args.output, comment_symbol="#")
