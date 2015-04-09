#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from os import path
from Tools.BLAST import BLASTPlus
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input fasta file")
parser.add_argument("-m", "--mask_file_name", action="store", dest="mask_filename",
                    help="Name of mask file")
parser.add_argument("-n", "--db_name", action="store", dest="db_name",
                    help="Name of blast db")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="format of input file")
args = parser.parse_args()

mask_filename = args.mask_filename if args.mask_filename else path.basename(args.input_file) + ".asnb"

BLASTPlus.make_blast_plus_db(args.input_file, mask_filename, args.db_name, format=args.format)

