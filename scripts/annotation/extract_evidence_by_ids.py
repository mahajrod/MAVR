#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import AUGUSTUS


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with AUGUSTUS evidence")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-d", "--id_file", action="store", dest="id_file", required=True,
                    help="File with ids to extract")
parser.add_argument("-m", "--mode", action="store", dest="mode", default="transcript",
                    help="Prefix of output files. Default - transcript")
args = parser.parse_args()

AUGUSTUS.extract_evidence_by_ids(args.input, args.id_file, args.output, mode=args.mode)
