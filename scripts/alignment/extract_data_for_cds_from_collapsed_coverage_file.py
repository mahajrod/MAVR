#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Bedtools import GenomeCov

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input collapsed coverage file")
parser.add_argument("-c", "--cds_bed", action="store", dest="cds_bed", required=True,
                    help="BED file with CDS coordinates")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output collapsed coverage file for CDS")
parser.add_argument("-s", "--skip_transcript_with_no_cds", action="store_true", dest="skip_transcript_with_no_cds",
                    help="Skip transcripts with no CDS. Default - false")

args = parser.parse_args()

GenomeCov.extract_data_for_cds_from_collapsed_coverage_file(args.input, args.cds_bed, args.output,
                                                            skip_transcript_with_no_cds=args.skip_transcript_with_no_cds)
