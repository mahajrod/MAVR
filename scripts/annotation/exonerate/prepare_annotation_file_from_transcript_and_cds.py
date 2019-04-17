#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import Exonerate


parser = argparse.ArgumentParser()

parser.add_argument("-t", "--transcript_file", action="store", dest="transcript_file", required=True,
                    help="File with sequences of transcripts")
parser.add_argument("-c", "--cds_file", action="store", dest="cds_file", required=True,
                    help="File with sequences of CDS")
parser.add_argument("-r", "--correspondence_file", action="store", dest="correspondence_file", required=True,
                    help="File with correspondence of transcripts to CDS")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input files. Allowed: fasta(default), genbank")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Verbose mode")

args = parser.parse_args()

Exonerate.prepare_annotation_file_from_transcript_and_cds(args.transcript_file, args.cds_file, args.correspondence_file,
                                                          args.output_prefix, format=args.format,
                                                          correspondence_key_column=0, correspondence_value_column=1,
                                                          verbose=args.verbose)
