#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from Tools.Annotation import AUGUSTUS

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fasta file with sequences")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output file")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use")

parser.add_argument("-s", "--species", action="store", dest="species", required=True,
                    help="Species to use as model")
parser.add_argument("-r", "--strand", action="store", dest="strand", default="both",
                    help="Strand to consider. Possible variants: both, forward, backward."
                         "Default: both")
parser.add_argument("-g", "--gene_model", action="store", dest="gene_model", default="complete",
                    help="Gene model to use. Possible variants:"
                         "partial      : allow prediction of incomplete genes at the sequence boundaries (default)"
                         "intronless   : only predict single-exon genes like in prokaryotes and some eukaryotes"
                         "complete     : only predict complete genes"
                         "atleastone   : predict at least one complete gene"
                         "exactlyone   : predict exactly one complete gene"
                         "Default: complete")

parser.add_argument("-e", "--other_options", action="store", dest="other_options",
                    help="Other augustus options")
parser.add_argument("-c", "--augustus_config_dir", action="store", dest="config_dir",
                    help="Augustus config dir")



args = parser.parse_args()

AUGUSTUS.threads = args.threads

AUGUSTUS.parallel_predict(args.species, args.input, args.output, strand=args.strand, gene_model=args.gene_model,
                          output_gff3=True, other_options=args.other_options, config_dir=args.config_dir)