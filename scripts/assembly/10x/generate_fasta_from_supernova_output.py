#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
import matplotlib

from RouToolPa.Tools.Assemblers import Supernova

parser = argparse.ArgumentParser()

parser.add_argument("-a", "--assembly_dir", action="store", dest="assembly_dir", required=True,
                    help="Directory with Supernova output, i.e /path/to/outs/assembly")
parser.add_argument("-m", "--min_length", action="store", dest="min_length", default=1000, type=int,
                    help="Minimum length of scaffold to output. Default: 1000")
parser.add_argument("-e", "--header_style", action="store", dest="header_style", default="full",
                    help="Fasta header style to use. Allowed: full(default), short")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output prefix")

args = parser.parse_args()

Supernova.generate_fasta(args.assembly_dir, args.out_prefix, min_length=args.min_length, header_style=args.header_style)
