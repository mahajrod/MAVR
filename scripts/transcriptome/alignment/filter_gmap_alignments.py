#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Alignment import GMAP

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input gff")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-c", "--min_coverage", action="store", dest="min_coverage", type=float,
                    help="Minimum coverage to retain transcripts. Default: not set")
parser.add_argument("-d", "--min_identity", action="store", dest="min_identity", type=float,
                    help="Minimum identity to retain transcripts. Default: not set")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Verbose output")
args = parser.parse_args()

GMAP.filter_gff(args.input, args.output, min_coverage=args.min_coverage,
                min_identity=args.min_identity, gene_feature="gene",
                transcript_feature="mRNA", exon_feature="exon", cds_feature="CDS", verbose=args.verbose)
