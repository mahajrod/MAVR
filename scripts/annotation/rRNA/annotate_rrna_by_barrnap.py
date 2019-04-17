#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Tools.Annotation import Barrnap

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_fasta", action="store", dest="input", required=True,
                    help="Input .fasta file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("-k", "--kingdom", action="store", dest="kingdom", default="euk",
                    help="Kingdom to use for rrna annotation. Allowed: bac, euk, mito, arc."
                         "Default: euk")
parser.add_argument("-l", "--length_cutoff", action="store", dest="length_cutoff", type=float,
                    default=0.8, help="Annotated rRNAs with length below cutoff will be treated as partial. "
                                      "Default: 0.8")
parser.add_argument("-r", "--reject_cutoff", action="store", dest="reject_cutoff", type=float,
                    default=0.001, help="Annotated rRNAs with length below cutoff will be rejected."
                                        "Default: 0.001")
parser.add_argument("-e", "--evalue_cutoff", action="store", dest="evalue_cutoff", type=float,
                    default=0.000001,
                    help="Hits with evalue below cutoff will be ignored. Default: 0.000001")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use. Default: 1")
parser.add_argument("-b", "--barrnap_dir", action="store", dest="barrnap_dir", default="",
                    help="Path to Barrnap directory")

args = parser.parse_args()

Barrnap.threads = args.threads
Barrnap.path = args.barrnap_dir
Barrnap.annotate_rrna(args.input, args.kingdom, args.output_prefix, length_cutoff=args.length_cutoff,
                      reject_cutoff=args.reject_cutoff, evalue_cutoff=args.evalue_cutoff)

