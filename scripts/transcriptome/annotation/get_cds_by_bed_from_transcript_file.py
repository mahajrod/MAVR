#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import SequenceRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with transcripts")
parser.add_argument("-o", "--output_file", action="store", dest="output", required=True,
                    help="Output file with CDS")
parser.add_argument("-c", "--coordinates_bed", action="store", dest="coordinates", required=True,
                    help="BED file with coordinates")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of transcript file. Default - fasta")
parser.add_argument("-p", "--parsing_mode", action="store", dest="parsing_mode", default="parse",
                    help="Parsing mode. Default - parse")
parser.add_argument("-z", "--zerobased", action="store_true", dest="zerobased", default=False,
                    help="Coordinates are zerobased. Default - False")

args = parser.parse_args()

SequenceRoutines.get_cds_by_bed_from_transcript_file(args.input, args.coordinates, args.output, zero_based=False,
                                                     transcript_format=args.format, parsing_mode="parse")
